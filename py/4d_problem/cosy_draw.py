#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

import sys, math
import os, shutil
#import commands
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append, zeros

fNom = array([870.9084099103394, 0.01625617474106655, 0.02579740767554017, 0.0020116413429212])
fNom = array([256.5588046909630, 0.2579740767554017E-001])
fNom = array([374.8469454403801, 0.2258636266940250])
#fNom = array([58.05270824959980, 0.7146266618279268E-001])
qNom = zeros(11)+1
qNew = array([1.00311026,1.02485471,0.86210439,0.71888156,1.06544479,0.98006624,0.70365153])
#qNew = array([0.77403563,0.90112332,1.24112697,1.95036434,1.80219943,0.67776703,0.84913223])
#qNew = array([-0.319612, 0.198440, 0.221305,-0.179820, 0.090087, 0.194583,-0.040515])

qNew = array([1.973248348851267,3.435842142158499,1.863635914537548,0.263583540331633,2.316445182494209,1.262226405304527,0.493986935739694])
PYGMO_DIR = '../../'
FOX_DIR = PYGMO_DIR + 'fox/'
COSY_DIR = '/mnt/simulations/secarml/COSY/'

def write_fox(qs=qNom, name=None, directory=FOX_DIR):
    input_len = len(qs)
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
#    [q1s, q2s, q3s, q4s, q5s, q6s, q7s, q8s, q9s, q10s, q11s] = qs
    if name==None:
        rand = rng()
    else:
        rand = name
    cosy_input = FOX_DIR + '20Ne1.18-3.5umCFoil_draw4d.fox'
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
    cmd = 'cosy ' + cosyFilename
#    failure, output = commands.getstatusoutput(cmd)
    output = commands.run(['cosy',cosyFilename], capture_output=True)
    stripped = output.stdout.strip()
#    print(stripped.split())
    resol = (stripped.split())[-2:]
    print(stripped)
    print(resol)
    for i in range(len(resol)):
        try: 
            resol[i] = float(resol[i])
        except:
            if i == 0:
                resol[i] = float(1e-6)
            else:
                resol[i] = float(1e6)
        if i == 0:
            resol[i] = fNom[i]/resol[i]
        else:    
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



