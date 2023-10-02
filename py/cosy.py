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

# define the nominal qvalue array (array is sent to cosy as a power of 2, i.e. 0 => 2^0 = 1 * nominal value)
qNom = zeros(19)+1
qNew = array([0.30961,1.0464,1.71728,0.70975,0.41581,0.51202,0.26017,1.30771,1.55515,0.22818,0.39962,0.30686,0.31115,0.72737,0.54035,1.22764,0.80459,3.90926,3.46931])
qNew = array([0.283413424223053,1.058069943743000,1.467197412744679,0.584709930968694,0.274509693545302,0.473285293926786,0.254461623764208,0.402837250182905,0.252802563833367, 0.21680981887451, 0.35877567899768, 0.42672424886372, 0.52992533228021, 0.49236325403560, 0.38959174127019,0.251246952092431,3.987305095918524,3.969383625652894,0.367206626871934])
qNew = array([0.2651,0.99709,1.69439,0.72239,0.44316,0.53659,0.26321,0.6078,0.28921,0.9526,0.27656,0.35334,0.30253,1.06311,0.75223,0.72,3.07214,3.83866,3.94119])
qNew = array([0.29933,1.13503,1.36588,0.55368,0.26,0.73098,0.28063,0.32898,0.75059,0.61074,0.57724,0.25125,0.43583,0.69772,0.6391,0.32272,3.99534,3.83737,0.72064])
qNew = array([0.973223,1.034567,0.965811,0.909692,0.883331,0.937646,0.687216,0.824256,0.972822,0.775758,0.980855,0.941599,0.892135,1.112539,1.216678,1.038781,0.907588,0.725516,1.167915])
qNew = array([0.979312, 0.979985, 1.003724, 0.976891, 0.908027])
#qNew = array([-1.691493, 0.065433, 0.780122, -0.494614, -1.266011, -0.965742, -1.942487, 0.387044, 0.637050, -2.131748, -1.323307, -1.704343, -1.684295, -0.459238, -0.888028, 0.295885, -0.313672, 1.966895, 1.794648])
#qNew = power(2, qNew)
#qNew = array([0.29484,1.04943,0.25012,1.52012,0.60846,0.29305,3.99557,0.50338,0.25181,3.99858,0.75819,0.65293,0.28756,0.2302,0.34915,0.30065,0.31636,0.7858,0.35506])
# set working DIR for PYGMO, FOX, COSY
PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
COSY_DIR = PYGMO_DIR + 'COSY10.0/'
SCRATCH_DIR = configs['scratch_dir']

on_fireside = secar_utils.check_fireside()
#if on_fireside:
#    FOX_DIR = SCRATCH_DIR + 'fox/'
#    isExist = os.path.exists(FOX_DIR)
#    if not isExist:
#       os.makedirs(FOX_DIR)

def cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2):
    # remove old cosy fox and lis file
    if os.path.exists(cosyFilename):
        os.remove(cosyFilename)
    if os.path.exists(lisFilename):
        os.remove(lisFilename)
    if os.path.exists(cosyFilename2):
        os.remove(cosyFilename2)
    if os.path.exists(lisFilename2):
        os.remove(lisFilename2)
    return 

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
    cosy_input = directory + fox_file
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
def cosyrun(qs, fNom, fox_name="SECAR_pg_Optics", directory_name=FOX_DIR):
    
    magnet_names = array(["Q1", "Q2", "B1", "B2", "HEX1", "Q3", "Q4", "Q5", "B3", "B4", "HEX2", "Q6", "Q7", "HEX3", "OCT1", "Q8", "Q9", "B5", "B6", "Q10", "Q11", "Q12", "Q13", "B7", "B8", "Q14", "Q15", "UMCP", "Viewer","DSSD"])
    # define the dimensions of the magnet chambers
    magnet_dims = array([[90,80],[140,102],[240,60],[240,60],[240,142],[220,142],[146,126],[102,102],[156,104],[240,104],[240,102],[280,110],[280,110],[165,115],[102,102],[100,100],[120,90],[148,66],[148,66],[180,96],[240,91],[140,140],[100,100],[130,60],[130,60],[100,100],[100,100],[75/1.41,75/1.41],[70,70],[64,64]])
    
    last_element = 'UMCP'
    last_element_index = np.where(magnet_names == last_element)[0][0]
#    print(last_element_index)

    # make fox file and get name
    fox_file = "{}.fox".format(fox_name)
    fox_DE_file = "{}_DE.fox".format(fox_name)
    cosyFilename, lisFilename = write_fox(qs,fox_file=fox_file, directory=directory_name)
    cosyFilename2, lisFilename2 = write_fox(qs,fox_file=fox_DE_file, directory=directory_name)
    
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
    n_objs = configs['n_obj'] + configs['n_con']
    objectives_constraints = []
    objectives_constraints = configs['objectives'] + configs['constraints']
    resol = zeros(n_objs) 

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
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
                return resol        
            fp1xdim_bool = False
        if fp2xdim_bool:
            try:
                fp2xdim = (float(split[i])*1000)
            except:
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
                return resol        
            fp2xdim_bool = False
        if fp3xdim_bool:
            try:
                fp3xdim = (float(split[i])*1000)
            except:
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
                return resol        
            fp3xdim_bool = False
        if beamspotsize_bool:
            try:
                beamspotsize = power(float(split[i]),0.5)
            except:
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
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
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
                return resol        
            fp2res_bool = False
        if fp3res_bool:
            try:
                fp3res = (float(split2[i])*1000)
            except:
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
                return resol        
            fp3res_bool = False
        if fp1res_bool:
            try:
                fp1res = (float(split2[i])*1000)
            except:
                resol = zeros(n_objs)+scale         
                cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
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

    # if no x-ydim values, just return outside constraints (all scale)
    if len(xdim) < len(magnet_dims) or len(ydim) < len(magnet_dims):
        resol = zeros(n_objs)+scale         
        cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
        return resol     
    dssd_x = abs(xdim[len(xdim)-1] * 1000)
    dssd_y = abs(ydim[len(ydim)-1] * 1000)
    beamspotsize = min(power(max(dssd_x/magnet_dims[np.where(magnet_names=='DSSD')][0][0], dssd_y/magnet_dims[np.where(magnet_names=='DSSD')][0][1],1),4), scale*10)
#    print(dssd_x, dssd_y, beamspotsize)
    if beamspotsize > scale:
        resol = zeros(n_objs)+scale         
        cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
        return resol 

    # cut off the dimensions at the UMCP
    xdim = xdim[:last_element_index+1]
    ydim = ydim[:last_element_index+1]
    magnet_dims = magnet_dims[:last_element_index+1]

    max_width = 0
    dict_resol = {}

    for i in range(len(magnet_dims)):
        if fp2xdim == 0 or fp3xdim == 0 or isnan(fp2xdim) or isnan(fp3xdim) or xdim2[i]==0 or ydim2[i]==0:
            resol = zeros(n_objs)+scale         
            cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
            return resol            
        # find xbound, ybound
        xbound = (abs(xdim[i])*1000)
        ybound = (abs(ydim[i])*1000)
        # if ratio of *bound to magnet_dim is larger than max_width, set max_width
        max_width = max(xbound/magnet_dims[i][0],ybound/magnet_dims[i][1],max_width)
        # if *bound more than 4x magnet_dim, or is nan, return outside constraints
#        print(magnet_names[i], xbound, ybound, magnet_dims[i][0],magnet_dims[i][1],max_width)
        if xbound > magnet_dims[i][0] * scale or ybound > magnet_dims[i][1] * scale or isnan(xbound) or isnan(ybound):
            resol = zeros(n_objs)+scale         
            cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
            return resol

#    max_width = min(power(max_width, 4), scale*10)
#    if max_width > scale*0.01:
#        resol = zeros(n_objs)+scale         
#        cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
#        return resol 
#    elif max_width > power(4.4,4):
#        max_width *= 100

    # if within constraints, set resol temporarily
    try:
        for i in configs['fixed_cosy_objs']:
            if i == 'FP1_res':
                dict_resol[i] = fp1xdim/fp1res
            if i == 'FP2_res':
                dict_resol[i] = fp2xdim/fp2res
            if i == 'FP3_res':
                dict_resol[i] = fp3xdim/fp3res
            if i == 'MaxBeamWidth':
                dict_resol[i] = max_width 
            if i == 'FP4_BeamSpot':
                dict_resol[i] = beamspotsize 
    except:
        resol = zeros(n_objs)+scale         
        cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
        return resol 

    for i in range(len(objectives_constraints)):
        if dict_resol[objectives_constraints[i]] < scale and dict_resol[objectives_constraints[i]] > 0:
#            print(objectives_constraints[i],dict_resol[objectives_constraints[i]])
            resol[i] = dict_resol[objectives_constraints[i]] / fNom[objectives_constraints[i]]
        else:
            resol = zeros(n_objs)+scale         
            cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)
            return resol 

    cleanup_fox(cosyFilename, lisFilename, cosyFilename2, lisFilename2)

    # return the objective values
    return (resol)
    

if __name__ == "__main__":
    # if running from console, just run the nominal setting
    fNom = configs['fNominal']
    print(fNom)
    print(cosyrun(qNew,fNom,configs['fox_name']))
    print(cosyrun(qNom,fNom,configs['fox_name']))

    fNom_keys = ["FP1_res","FP2_res","FP3_res","MaxBeamWidth","FP4_BeamSpot"]
    fNom_values = array([2.3849360856494263, .10961548781662407, .5108029152516118, 1.6251646888022029, 1.0])
    fNom = {}
    for i in range(len(fNom_keys)):
        fNom[fNom_keys[i]] = fNom_values[i]
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


