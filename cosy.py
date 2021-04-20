#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python3.9

import sys, math
import os, shutil
#import commands
import subprocess as commands
import re
from numpy.random import random as rng
from numpy import array, append

qNom = array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )
fNom = array([-870.9084099103394, 0.01625617474106655, 0.02579740767554017, 0.0020116413429212])
qNom = array([1,1,1,1,1,1,1] )
qNew = array([1.00311026,1.02485471,0.86210439,0.71888156,1.06544479,0.98006624,0.70365153])
#qNew = array([0.77403563,0.90112332,1.24112697,1.95036434,1.80219943,0.67776703,0.84913223])
#qNew = array([-0.319612, 0.198440, 0.221305,-0.179820, 0.090087, 0.194583,-0.040515])

# Function that runs cosy given field gradients and outputs resolution at FP3. 
# Output is written in file temp-results
def cosyrun(qs=qNom):
    input_len = len(qs)
    if (len(qNom)-input_len>0):
        for i in range(len(qNom)-input_len):
            qs = append(qs,qNom[i+input_len])
    [q1s, q2s, q3s, q4s, q5s, q6s, q7s] = qs
    rand = rng()
    cosy_input = '20Ne1.18-3.5umCFoil.fox'
    text = None
    with open(cosy_input, 'r') as f:
        text = f.readlines()

    start_line = 563
    for i in range(len(qs)):
        text[i+start_line] = "Q{0}:= {1};\n".format(i+1,qs[i])

    cosyFilename = 'pygmoCosy'+str(rand)+'.fox'
    while os.path.exists(cosyFilename):
        rand = rng()
        cosyFilename = 'pygmoCosy'+str(rand)+'.fox'
    lisFilename = 'pygmoCosy'+str(rand)+'.lis'
    # creating input file
    with open(cosyFilename, 'w') as f:
        f.writelines(text)
    
    #Removing files from older runs
#    failure, output = commands.getstatusoutput(cmd)
    
    #Run file
    cmd = 'cosy ' + cosyFilename
#    failure, output = commands.getstatusoutput(cmd)
    output = commands.run(['cosy',cosyFilename], capture_output=True)
    stripped = output.stdout.strip()
#    print(stripped.split())
    resol = (stripped.split())
    print(resol)
    for i in range(len(resol)):
        resol[i] = float(resol[i])
        if i == 0:
            if resol[i] < 0:
                resol[i] = fNom[i]/resol[i]
            else:
                resol[i] = max(resol[i],100000)
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
    print(cosyrun(qNew))



