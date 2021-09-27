#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

import sys, math
import os, shutil
#import commands
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append, zeros, power

# Q1, Q2, B1, B2, HEX1, Q3, Q4, Q5, B3, B4, HEX2, Q6, Q7, HEX3, OCT1, Q8, Q9, B5, B6, Q10, Q11
magnet_dims = array([[90,80],[162,102],[240,60],[240,60],[242,142],[242,142],[146,126],[102,102],[156,104],[156,104],[240,102],[290,110],[290,110],[165,115],[102,102],[121,101],[156,90],[148,66],[148,66],[222,96],[282,91]])
fNom = array([870.9084099103394, 0.01625617474106655, 0.02579740767554017, 0.0020116413429212])
fNom = array([256.5588046909630, 0.2579740767554017E-001])
fNom = array([255.8823693915960, 0.2579740767554017E-001,374.8469454403801, 0.2258636266940250])
#fNom = array([58.05270824959980, 0.7146266618279268E-001])
qNom = zeros(11)+1
qNew = array([1.00311026,1.02485471,0.86210439,0.71888156,1.06544479,0.98006624,0.70365153])
#qNew = array([0.77403563,0.90112332,1.24112697,1.95036434,1.80219943,0.67776703,0.84913223])
#qNew = array([-0.319612, 0.198440, 0.221305,-0.179820, 0.090087, 0.194583,-0.040515])
qNew = array([-0.281811,-0.155926,-0.052414,-0.157868,-0.007227, 0.047082, 0.119611,-0.211220,-0.040499, 0.126084, 0.020534])
qNew = power(zeros(11)+2,qNew)
PYGMO_DIR = '../../'
FOX_DIR = PYGMO_DIR + 'fox/'
COSY_DIR = '/mnt/simulations/secarml/COSY/'

def write_fox(qs=qNom, name=None, directory=FOX_DIR):
    input_len = len(qs)
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    [q1s, q2s, q3s, q4s, q5s, q6s, q7s, q8s, q9s, q10s, q11s] = qs
    if name==None:
        rand = rng()
    else:
        rand = name
    cosy_input = FOX_DIR + '20Ne1.18-3.5umCFoil_ydims.fox'
    text = None
    with open(cosy_input, 'r') as f:
        text = f.readlines()

    start_line = 0 
    for i in range(len(text)):
        if "SET MAGNETS" in text[i]:
            start_line = i+1
            print(text[i], text[start_line])
    for i in range(len(qs)):
        text[i+start_line] = "Q{0}:= {1};\n".format(i+1,qs[i])

    cosyFilename = directory + 'pygmoCosy'+str(rand)+'.fox'
#    print(cosyFilename)
    while os.path.exists(cosyFilename):
        rand = rng()
        cosyFilename = directory + 'pygmoCosy'+str(rand)+'.fox'
    lisFilename = directory + 'pygmoCosy'+str(rand)+'.lis'
    # creating input file
    with open(cosyFilename, 'w') as f:
        f.writelines(text)
    return cosyFilename, lisFilename


# Function that runs cosy given field gradients and outputs resolution at FP3. 
# Output is written in file temp-results
def cosyrun(qs=qNom):

    cosyFilename, lisFilename = write_fox(qs)
    
    #Removing files from older runs
#    failure, output = commands.getstatusoutput(cmd)
    
    #Run file
    cmd = COSY_DIR + 'cosy'
#    failure, output = commands.getstatusoutput(cmd)
    output = commands.run([cmd ,cosyFilename], capture_output=True)
    stripped = output.stdout.strip().decode('utf8','strict')
    xmax, xmin = [], []
    ymax, ymin = [], []
    fp2res, fp2espread = 0, 0
    fp3res, fp3espread = 0, 0
    xmax_bool, xmin_bool = False, False
    ymax_bool, ymin_bool = False, False
    fp2res_bool, fp2espread_bool = False, False
    fp3res_bool, fp3espread_bool = False, False
    split = stripped.split()
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
    resol1 = (stripped.split())[:2]
    resol2 = (stripped.split())[-2:]
#    resol = [resol1[0],resol1[1],resol2[0],resol2[1]]
    resol = [fp2res,fp2espread,fp3res,fp3espread]
    print(split)
    print(xmax, ymax)
    print(xmin, ymin)
    print(resol)
    for i in range(len(resol)):
        if resol[i] > 0: 
            resol[i] = float(resol[i])
        else:
            if i == 0 or i == 2:
                resol[i] = float(1e-9)
            else:
                resol[i] = float(1e9)
        if i == 0 or i == 2:
            resol[i] = fNom[i]/resol[i]
        else:    
            resol[i] = resol[i]/fNom[i]
    for i in range(len(magnet_dims)):
        if len(xmax) == 0 or len(xmin) == 0 or len(ymax) == 0 or len(ymin) == 0:
            resol = array([1e9,1e9,1e9,1e9])         
            break            
        xbound = max(abs(xmax[i] * 1000),abs(xmin[i] * 1000))
        ybound = max(abs(ymax[i] * 1000),abs(ymin[i] * 1000))
        if xbound > magnet_dims[i][0] or ybound > magnet_dims[i][1]:
            resol = array([1e9,1e9,1e9,1e9])         
            break
    print(resol)            
    commands.run(['rm','-f',cosyFilename])
    commands.run(['rm','-f',lisFilename])
#    try:
#        f = open(tempresFilename,'r')  # Opening temp file with results
#        resol = f.readline()
#        f.close()
#    except (OSError, IOError) as e:
#        print('q1 = %.6f q2 = %.6f q3 = %.6f q4 = %.6f q5 = %.6f q6 = %.6f q7 = %.6f'  %(q1s, q2s, q3s, q4s, q5s, q6s, q7s) )   
#        resol = 0
#    f = open('results.txt','a')  # Writing results file: magnet field, resolution
#    f.write( '{0:d} {1:.6f} {2:.6f} {3:.6f} {4:.6f} {5:.6f} {6:.6f} {7:.6f} {8:.1f}\n' .format(b, q1s, q2s, q3s, q4s, q5s, q6s, q7s, float(resol)) )
#    f.close()

#    print('Finished bee %d'%b)
    return (resol)
    

if __name__ == "__main__":
#    print(cosyrun(qNom))
    print(cosyrun(qNew))



