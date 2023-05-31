#!/usr/bin/env python3

#import commands
import sys, math
import os, shutil
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append, zeros, power, isnan, divide
import numpy as np
import pandas as pd
import timeit

import secar_utils as secar_utils
configs = secar_utils.load_configs()

magnet_names = array(["Q1", "Q2", "B1", "B2", "HEX1", "Q3", "Q4", "Q5", "B3", "B4", "HEX2", "Q6", "Q7", "HEX3", "OCT1", "Q8", "Q9", "B5", "B6", "Q10", "Q11", "Q12", "Q13", "B7", "B8", "Q14", "Q15", "UMCP", "Viewer"])
# define the dimensions of the magnet chambers
magnet_dims = array([[90,80],[140,102],[240,60],[240,60],[240,142],[220,142],[146,126],[102,102],[156,104],[156,104],[240,102],[280,110],[280,110],[165,115],[102,102],[100,100],[120,90],[148,66],[148,66],[180,96],[240,91],[140,140],[100,100],[130,60],[130,60],[100,100],[100,100],[75/1.41,75/1.41],[70,70]])

# define the nominal qvalue array (array is sent to cosy as a power of 2, i.e. 0 => 2^0 = 1 * nominal value)
qNom = zeros(19)+1
qNew = array([0.9995845752046637,
1.2220848121694863,
0.7288452110421395,
0.5452658512360776,
0.25893655117638975,
0.8480128005770239,
0.5387917381151589,
0.9130364420922832,
0.987050172031228,
0.41523318638126994,
0.9503558977890336,
0.9020902091901309,
0.7972568272495952,
0.9012788345387545,
0.6349413818999178,
0.8762447428324577,
0.3558285257840111,
0.4032267640375105,
0.2894016712233756])

# set working DIR for PYGMO, FOX, COSY
PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
COSY_DIR = PYGMO_DIR + 'COSY10.0/'
SCRATCH_DIR = configs['scratch_dir']

on_fireside = secar_utils.check_fireside()
if on_fireside:
    FOX_DIR = SCRATCH_DIR + 'fox/'
    isExist = os.path.exists(FOX_DIR)
    if not isExist:
       os.makedirs(FOX_DIR)

# write the qvalue array to a fox file for cosy to run
def write_fox(qs=qNom, name=None, directory=FOX_DIR, fox_file='SEC_neutrons_WF_14m_v1.fox'):
    # how many magnets to set
    input_len = len(qs)
    # if less than qNom, use qNom values (0) for remainder
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    # get rand number for the temporary file
    if name==None:
        rand = rng()
    else:
        rand = name
    # get the fox file for the simulation and change the qvalues
    cosy_input = FOX_DIR + fox_file
    text = None
    with open(cosy_input, 'r') as f:
        text = f.readlines()

    start_line = 0 
    # find where in the file to change
    for i in range(len(text)):
        # dummy text line in the fox file that the code searches for to find the qvalue setting
        if "SET MAGNETS" in text[i]:
            start_line = i+1
#            print(text[i], text[start_line])

    # change the q setttings
    magnet_i = 0
    for i in range(len(qs)):
        magnet_i += 1
        if i < 15:
            text[i+start_line] = "Q{0}_SC:= {1};\n".format(magnet_i,qs[i])
        elif i < 18: 
            magnet_i = i - 14
            text[i+start_line] = "H{0}_SC:= {1};\n".format(magnet_i,qs[i])
        else: 
            magnet_i = 1
            text[i+start_line] = "O{0}_SC:= {1};\n".format(magnet_i,qs[i])

    # temporary output file
    cosyFilename = directory + 'pygmoCosy'+str(rand)+'.fox'
    # if by some ridiculous chance we pick the same number as another iteration, find new number
    while os.path.exists(cosyFilename):
        rand = rng()
        cosyFilename = directory + 'pygmoCosy'+str(rand)+'.fox'
    lisFilename = directory + 'pygmoCosy'+str(rand)+'.lis'

    # write out the fox file
    with open(cosyFilename, 'w') as f:
        f.writelines(text)
    # return filenames for cosy to find
    return cosyFilename, lisFilename


# Function that runs cosy 
# check problem.py and optimize.py for more on how the operation works
# idea is that evolve() in optimize.py will call fitness() in problem.py
#   fitness() will pass a set of magnet factors (between 0.5 and 2.0) to cosyrun
#   and will return the resulting objectives to the algorithm
def cosyrun(qs, fNom, fox_name="SECAR_pg_Optics"):

    # make fox file and get name
    fox_file = "{}.fox".format(fox_name)
    fox_DE_file = "{}_DE.fox".format(fox_name)
    cosyFilename, lisFilename = write_fox(qs,fox_file=fox_file)
    cosyFilename2, lisFilename2 = write_fox(qs,fox_file=fox_DE_file)
    
    #Run cmd
    cmd = COSY_DIR + 'cosy'
    # timer for diagnostics
    startTime = timeit.default_timer()
    # run cosy 
    output = commands.run([cmd ,cosyFilename], capture_output=True)
    # print time
#    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))
    # timer for diagnostics
    startTime = timeit.default_timer()
    # run cosy2 
    output2 = commands.run([cmd ,cosyFilename2], capture_output=True)
    # print time
#    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

    # get output and now convert into the necessary values to return to pygmo
    stripped = output.stdout.strip().decode('utf8','strict')
    split = stripped.split()
#    print(split)
    # get output and now convert into the necessary values to return to pygmo
    stripped2 = output2.stdout.strip().decode('utf8','strict')
    split2 = stripped2.split()
#    print(split2)

    # initiate all variables to be read, and bools for the reader to check
    xdim, ydim = [], []
    xdim2, ydim2 = [], []
    fp1res, fp1espread, fp1xdim = 0, 0, 0
    fp2res, fp2espread, fp2xdim = 0, 0, 0
    fp3res, fp3espread, fp3xdim = 0, 0, 0
    beamspotsize = 0
    xdim_bool, ydim_bool = False, False
    fp1res_bool, fp1espread_bool, fp1xdim_bool = False, False, False
    fp2res_bool, fp2espread_bool, fp2xdim_bool = False, False, False
    fp3res_bool, fp3espread_bool, fp3xdim_bool = False, False, False
    beamspotsize_bool = False

    # scale factor to account for the beam spot issue
    #   even the nominal setting is outside the bounds...
    scale = 1e9 
    # setup value to be returned, here 4 different objectives
    #objs = configs['n_obj']
    objs = 5
    resol = zeros(objs) 

    # my method could probably be better optimized but this works and is mostly straight-forward
    #   idea is just that keywords are output for each variable, e.g. "Xdim" for the Xdim for a magnet dimension
    #   the below code will check if the bool is true and if so write the corresponding value and false the bool
    #   in total we parse through all the cosy output and find all the different variables
    for i in range(len(split)):
        if xdim_bool:
            xdim.append(float(split[i]))
            xdim_bool = False
        if ydim_bool:
            ydim.append(float(split[i]))
            ydim_bool = False
        if fp1xdim_bool:
            try:
                fp1xdim = (float(split[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp1xdim_bool = False
        if fp2xdim_bool:
            try:
                fp2xdim = (float(split[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp2xdim_bool = False
        if fp3xdim_bool:
            try:
                fp3xdim = (float(split[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp3xdim_bool = False
        if beamspotsize_bool:
            try:
                beamspotsize = power(float(split[i]),0.5)
            except:
                resol = zeros(objs)+scale         
                return resol        
            beamspotsize_bool = False
        if split[i].strip() == "Xdim":
            xdim_bool = True
        if split[i].strip() == "Ydim":
            ydim_bool = True
        if split[i].strip() == "FP1Xdim":
            fp1xdim_bool = True
        if split[i].strip() == "FP2Xdim":
            fp2xdim_bool = True
        if split[i].strip() == "FP3Xdim":
            fp3xdim_bool = True
        if split[i].strip() == "BeamSpotSize":
            beamspotsize_bool = True
    for i in range(len(split2)):
        if fp2res_bool:
            try:
                fp2res = (float(split2[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp2res_bool = False
        if fp3res_bool:
            try:
                fp3res = (float(split2[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp3res_bool = False
        if fp1res_bool:
            try:
                fp1res = (float(split2[i])*1000)
            except:
                resol = zeros(objs)+scale         
                return resol        
            fp1res_bool = False
        if xdim_bool:
            xdim2.append(float(split2[i]))
            xdim_bool = False
        if ydim_bool:
            ydim2.append(float(split2[i]))
            ydim_bool = False
        if split2[i].strip() == "FP1DE":
            fp1res_bool = True
        if split2[i].strip() == "FP2DE":
            fp2res_bool = True
        if split2[i].strip() == "FP3DE":
            fp3res_bool = True
        if split2[i].strip() == "Xdim":
            xdim_bool = True
        if split2[i].strip() == "Ydim":
            ydim_bool = True

    max_width = 0

    for i in range(len(magnet_dims)):
        # if no x-ydim values, just return outside constraints (all scale)
        if len(xdim) < len(magnet_dims) or len(ydim) < len(magnet_dims):
            resol = zeros(objs)+scale         
            return resol     
        if fp2xdim == 0 or fp3xdim == 0 or isnan(fp2xdim) or isnan(fp3xdim) or xdim2[i]==0 or ydim2[i]==0:
            resol = zeros(objs)+scale         
            return resol            
        # find xbound, ybound
        xbound = (abs(xdim[i])*1000)
        ybound = (abs(ydim[i])*1000)
        # if ratio of *bound to magnet_dim is larger than max_width, set max_width
        max_width = max(xbound/magnet_dims[i][0],ybound/magnet_dims[i][1],max_width)
        # if *bound more than 4x magnet_dim, or is nan, return outside constraints
        if xbound > magnet_dims[i][0] * scale or ybound > magnet_dims[i][1] * scale or isnan(xbound) or isnan(ybound):
            resol = zeros(objs)+scale         
            return resol
        # if within constraints, set resol temporarily
        try:
            resol = [fp1xdim/fp1res,fp2xdim/fp2res,fp3xdim/fp3res,max_width,beamspotsize]
        except:
            resol = zeros(objs)+scale         
            return resol 

    scale = 1e9
#    print(divide(resol, fNom))
#    print(divide(resol, fNom)[1:])
#    print(np.all(divide(resol, fNom)[1:] < scale))
    # if within constraints, set resol as a ratio to nominal
    if np.all(divide(resol, fNom)[1:] < scale):
        for i in range(len(resol)):
            # make sure we are working with positive numbers
            if resol[i] > 0: 
                resol[i] = float(resol[i])
            else:
                resol = zeros(objs)+scale         
                return resol
            # we want to minimize max_width and beamspotsize, so just take resol/fNom
            resol[i] = resol[i]/fNom[i]
    else:
        resol = zeros(objs)+scale         

    # remove old cosy fox and lis file
    os.remove(cosyFilename)
    os.remove(lisFilename)
    os.remove(cosyFilename2)
    os.remove(lisFilename2)

    # return the objective values
    return (resol)
    

if __name__ == "__main__":
    # if running from console, just run the nominal setting
    fNom = configs['fNominal']
    print(cosyrun(qNom,fNom,configs['fox_name']))
    print(cosyrun(qNew,fNom,configs['fox_name']))

    fNom = array([2.3849360856494263, .10961548781662407, .5108029152516118, 1.6251646888022029, 0.12574090408565933])
    print(cosyrun(qNom,fNom,"SECAR_an_Optics"))
    # if running from console, just run the assigned setting
#    df = pd.read_csv('results_280/magnet_factors.csv')
#    print(df)
#    for i in range(df.shape[0]):
#        print(cosyrun(array(df.iloc[i,:19])))
#    PROFILES_PATH = "./"
#    plot_i = 1
#    write_fox(qNew, str(plot_i), PROFILES_PATH , 'SECAR_an_Optics_draw.fox')
##    write_fox(qNew, str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE_draw.fox')
#    write_fox(qNew, str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE2.fox')


