#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python3.9

import sys, math
import os, shutil
#import commands
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append, zeros, isnan

fNom = array([256.5588046909630, 0.2579740767554017E-001])
fNom = array([255.8823693915960, 0.2579740767554017E-001,374.8469454403801, 0.2258636266940250])
#fNom = array([255.8823693915960,0.2579740767554017E-001,373.3145173506128,0.2258636266940250,386.4836031051735,0.3032627082174342])
qNom = zeros(15)+1

PYGMO_DIR = '/mnt/home/herman67/cosy/pygmo/'
FOX_DIR = PYGMO_DIR + 'fox/'

def write_fox(qs=qNom, name=None, directory=FOX_DIR):
    input_len = len(qs)
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    [q1s, q2s, q3s, q4s, q5s, q6s, q7s, q8s, q9s, q10s, q11s, q12s, q13s, q14s, q15s] = qs
    if name==None:
        rand = rng()
    else:
        rand = name
    cosy_input = FOX_DIR + '20Ne1.18-3.5umCFoil_4f_10WC.fox'
    text = None
    with open(cosy_input, 'r') as f:
        text = f.readlines()

    start_line = 561
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
    cmd = 'cosy ' + cosyFilename
#    failure, output = commands.getstatusoutput(cmd)
    output = commands.run(['cosy',cosyFilename], capture_output=True)
    stripped = output.stdout.strip()
#    print(stripped.split())
    resol1 = (stripped.split())[:2]
#    resol2 = (stripped.split())[-4:-2]
    resol3 = (stripped.split())[-2:]
    resol = [resol1[0],resol1[1],resol3[0],resol3[1]]
    print(stripped)
    print(resol)
    for i in range(len(resol)):
        try: 
            resol[i] = float(resol[i])
        except:
            if i in [0,2,4]:
                resol[i] = float(1e-9)
            else:
                resol[i] = float(1e9)
        if i in [0,2,4]:
            resol[i] = fNom[i]/resol[i]
            if isnan(resol[i]):
                resol[i] = float(1e-9)
                resol[i] = fNom[i]/resol[i]
        else:    
            resol[i] = resol[i]/fNom[i]
            if isnan(resol[i]):
                resol[i] = float(1e9)
                resol[i] = resol[i]/fNom[i]
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
    print(cosyrun(qNom))



