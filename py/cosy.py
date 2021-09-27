#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python


#import commands
import sys, math
import os, shutil
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append, zeros, power
import timeit

# Q1, Q2, B1, B2, HEX1, Q3, Q4, Q5, B3, B4, HEX2, Q6, Q7, HEX3, OCT1, Q8, Q9, B5, B6, Q10, Q11, Q12, Q13, B7, B8, Q14, Q15
# define the dimensions of the magnet chambers
magnet_dims = array([[90,80],[140,102],[240,60],[240,60],[240,142],[220,142],[146,126],[102,102],[156,104],[156,104],[240,102],[280,110],[280,110],[165,115],[102,102],[100,100],[120,90],[148,66],[148,66],[180,96],[240,91],[140,140],[100,100],[130,60],[130,60],[100,100],[100,100]])

# define the nominal values for the objective function
fNom = array([215.7997411346150, 237.3677022535934, 0.5090584934220124, 0.05151361676274798])

# define the nominal qvalue array (array is sent to cosy as a power of 2, i.e. 0 => 2^0 = 1 * nominal value)
qNom = zeros(11)+1

# define a non nominal qvalue array, if checking the values
qNew = array([1.029116370504024,1.027586221796089,0.948009341322751,0.846517309969750,0.888012998083179,0.989275012835886,0.919427494202570,0.977212813272792,0.979037350285672, 0.87447623778262, 1.00129050151537])

# set working DIR for PYGMO, FOX, COSY
PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
#hpcc servers
COSY_DIR = '/mnt/home/herman67/cosy/COSY10.0/'
#nscl servers
#COSY_DIR = '/mnt/simulations/secarml/COSY/'

# write the qvalue array to a fox file for cosy to run
def write_fox(qs=qNom, name=None, directory=FOX_DIR):
    # how many magnets to set
    input_len = len(qs)
    # if less than qNom, use qNom values (0) for remainder
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    [q1s, q2s, q3s, q4s, q5s, q6s, q7s, q8s, q9s, q10s, q11s] = qs
    # get rand number for the temporary file
    if name==None:
        rand = rng()
    else:
        rand = name
    # get the fox file for the simulation and change the qvalues
    cosy_input = FOX_DIR + '20Ne1.18-3.5umCFoil_ydims.fox'
    text = None
    with open(cosy_input, 'r') as f:
        text = f.readlines()

    start_line = 0 
    # find where in the file to change
    for i in range(len(text)):
        # dummy text line in the fox file that the code searches for to find the qvalue setting
        if "SET MAGNETS" in text[i]:
            start_line = i+1
            print(text[i], text[start_line])
    # change the q setttings
    for i in range(len(qs)):
        text[i+start_line] = "Q{0}:= {1};\n".format(i+1,qs[i])

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
def cosyrun(qs=qNom):

    # make fox file and get name
    cosyFilename, lisFilename = write_fox(qs)
    
    #Run cmd
    cmd = COSY_DIR + 'cosy'
    # timer for diagnostics
    startTime = timeit.default_timer()
    # run cosy 
    output = commands.run([cmd ,cosyFilename], capture_output=True)
    # print time
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

    # get output and now convert into the necessary values to return to pygmo
    stripped = output.stdout.strip().decode('utf8','strict')
    split = stripped.split()

    # initiate all variables to be read, and bools for the reader to check
    xmax, xmin = [], []
    ymax, ymin = [], []
    fp2res, fp2espread = 0, 0
    fp3res, fp3espread = 0, 0
    beamspotsize = 0
    xmax_bool, xmin_bool = False, False
    ymax_bool, ymin_bool = False, False
    fp2res_bool, fp2espread_bool = False, False
    fp3res_bool, fp3espread_bool = False, False
    beamspotsize_bool = False

    # my method could probably be better optimized but this works and is mostly straight-forward
    #   idea is just that keywords are output for each variable, e.g. "Xmax" for the Xmax for a magnet dimension
    #   the below code will check if the bool is true and if so write the corresponding value and false the bool
    #   in total we parse through all the cosy output and find all the different variables
    for i in range(len(split)):
        if xmax_bool:
            xmax.append(float(split[i]))
            xmax_bool = False
        if xmin_bool:
            xmin.append(float(split[i]))
            xmin_bool = False
        if ymax_bool:
            ymax.append(float(split[i]))
            ymax_bool = False
        if ymin_bool:
            ymin.append(float(split[i]))
            ymin_bool = False
        if fp2res_bool:
            fp2res = (float(split[i]))
            fp2res_bool = False
        if fp2espread_bool:
            fp2espread = (float(split[i]))
            fp2espread_bool = False
        if fp3res_bool:
            fp3res = (float(split[i]))
            fp3res_bool = False
        if fp3espread_bool:
            fp3espread = (float(split[i]))
            fp3espread_bool = False
        if beamspotsize_bool:
            beamspotsize = power(float(split[i]),0.5)
            beamspotsize_bool = False
        if split[i].strip() == "Xmax":
            xmax_bool = True
        if split[i].strip() == "Xmin":
            xmin_bool = True
        if split[i].strip() == "Ymax":
            ymax_bool = True
        if split[i].strip() == "Ymin":
            ymin_bool = True
        if split[i].strip() == "FP2Res":
            fp2res_bool = True
        if split[i].strip() == "FP2Espread":
            fp2espread_bool = True
        if split[i].strip() == "FP3Res":
            fp3res_bool = True
        if split[i].strip() == "FP3Espread":
            fp3espread_bool = True
        if split[i].strip() == "BeamSpotSize":
            beamspotsize_bool = True


    max_width = 0
    # setup value to be returned, here 4 different objectives
    resol = zeros(4) 
    for i in range(len(magnet_dims)):
        # if no x-ymax/min values, just return an error run (all 1e9)
        if len(xmax) == 0 or len(xmin) == 0 or len(ymax) == 0 or len(ymin) == 0:
            resol = array([1e9,1e9,1e9,1e9])         
            break            
        # find  
        xbound = max(abs(xmax[i] * 1000),abs(xmin[i] * 1000))
        ybound = max(abs(ymax[i] * 1000),abs(ymin[i] * 1000))
        max_width = max(xbound/magnet_dims[i][0],ybound/magnet_dims[i][1],max_width)
        print(i, xbound, magnet_dims[i][0], ybound, magnet_dims[i][1])
        if xbound > magnet_dims[i][0] or ybound > magnet_dims[i][1]:
            print(xbound, magnet_dims[i][0], ybound, magnet_dims[i][1])
            resol = array([1e9,1e9,1e9,1e9])         
            break
        resol = [fp2res,fp3res,max_width,beamspotsize]
    print(resol)
    if max(resol) < 1e9:
        for i in range(len(resol)):
            if resol[i] > 0: 
                resol[i] = float(resol[i])
            else:
                if i == 0 or i == 1:
                    resol[i] = float(1e-9)
                else:
                    resol[i] = float(1e9)
            if i == 0 or i == 1:
                resol[i] = fNom[i]/resol[i]
            else:    
                resol[i] = resol[i]/fNom[i]
    print(resol)            
    commands.run(['rm','-f',cosyFilename])
    commands.run(['rm','-f',lisFilename])

    return (resol)
    

if __name__ == "__main__":
    # if running from console, just run the nominal setting
    print(cosyrun(qNom))
    # if running from console, just run the assigned setting
#    print(cosyrun(qNew))


