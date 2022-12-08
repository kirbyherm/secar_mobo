#!/mnt/simulations/secarml/soft/anaconda3/bin/python
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

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
from numpy import array, append, zeros, power, isnan
import timeit
import pandas as pd


magnet_names = array(["Q1", "Q2", "B1", "B2", "HEX1", "Q3", "Q4", "Q5", "B3", "B4", "HEX2", "Q6", "Q7", "HEX3", "OCT1", "Q8", "Q9", "B5", "B6", "Q10", "Q11", "Q12", "Q13", "B7", "B8", "Q14", "Q15", "UMCP", "Viewer"])
# define the dimensions of the magnet chambers
magnet_dims = array([[90,80],[140,102],[240,60],[240,60],[240,142],[220,142],[146,126],[102,102],[156,104],[156,104],[240,102],[280,110],[280,110],[165,115],[102,102],[100,100],[120,90],[148,66],[148,66],[180,96],[240,91],[140,140],[100,100],[130,60],[130,60],[100,100],[100,100],[75/1.41,75/1.41],[70,70]])

# define the nominal values for the objective function
fNom = array([245.5333762546184, 256.5533534865096, 1.016965710603861, 0.0497233197451071])
fNom = array([0.02285401532682956, 0.04181594290692345, 3.422466427009127, 0.27344973981231574, 0.05])
fNom = array([1086.7911810119049, 1258.9642382235916, 1155.6819246495133, 3.422466427009127, 0.27344973981231574])
fNom = array([1431.8410759523508, 821.7565325150232, 650.6352599978524, 0.934467870935426, 0.03972091942829642])
fNom = array([2384.9360856494263, 1058.1013973315412, 1260.5797906816008, 4.129531004594597, 0.32301378801056696])
fNom = array([2384.9360856494263, 109.61548781662407, 510.8029152516118, 1.6251646888022029, 0.12574090408565933])
fNom = array([123.09714342469559, 0.14485759961115405, 0.1050839232170788, 1.017425650763311, 0.0502221220191954])
# define the nominal qvalue array (array is sent to cosy as a power of 2, i.e. 0 => 2^0 = 1 * nominal value)
qNom = zeros(19)+1

# define a non nominal qvalue array, if checking the values
qNew = array([1.0371301857113335,1.4897519431921593,0.5402003843384104,0.6080163749223835,0.5965351874518491,0.5279178522813484,0.8474952322221544,0.8290931192132953,0.7350223112146984,0.5049139345530922,0.969681779928563 ,0.8465270119223961,0.7261232553654523,0.6805787940919176,0.6772214286022437,1.6737045402403927,1.3151418622198896,0.8914897696929639,0.6144362243855045])
#qNew = array([0.9231306600055562, 1.3296051142151415, 2.0103962914793567, 1.2425012100959048, 1.3524328502381404, 0.399121814584066, 0.23041908744123996, 0.1490014243074733, 0.19435642454669186, 2.93870928376383, 1.0713201915318586, 0.5435953805162628, 0.3525535276959007, 1.1428457798432394, 0.9097919998417044, 2.6821630795610516, 3.7961922235155785, 0.944569380681556, 0.2538890733802974])
qNew = array([1.0370,1.3988,0.5838,0.6992,0.6577,0.6059,0.6713,0.6439,0.5136,0.5630,0.8576,0.7022,0.5073,1.0006,0.7944,0.5916,1.7621,0.6079,0.5481])
# set working DIR for PYGMO, FOX, COSY
qNew = power(zeros(19)+2,array([-4.718300e-01, 2.197192e-01,-7.287368e-01,-1.685577e+00,-1.434873e+00,-2.522156e-01,-2.970157e-01, 1.631719e-01, 8.892526e-01,-2.278044e+00,-1.506612e+00,-7.736163e-01,-9.422170e-01,-1.111204e+00, 1.062845e+00,-1.917745e+00, 7.777601e-01,-1.195297e+00,-9.878368e-01]))

PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
#hpcc servers
#COSY_DIR = '/mnt/home/herman67/cosy/COSY10.0/'
#nscl servers
COSY_DIR = PYGMO_DIR + 'COSY10.0/'

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
def cosyrun(qs=qNom):

    # make fox file and get name
    cosyFilename, lisFilename = write_fox(qs,fox_file="SECAR_pg_Optics.fox")
    cosyFilename2, lisFilename2 = write_fox(qs,fox_file="SECAR_pg_Optics_DE.fox")
    
    #Run cmd
    cmd = COSY_DIR + 'cosy'
    # timer for diagnostics
    startTime = timeit.default_timer()
    # run cosy 
    output = commands.run([cmd ,cosyFilename], capture_output=True)
    # print time
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))
    # timer for diagnostics
    startTime = timeit.default_timer()
    # run cosy2 
    output2 = commands.run([cmd ,cosyFilename2], capture_output=True)
    # print time
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

    # get output and now convert into the necessary values to return to pygmo
    stripped = output.stdout.strip().decode('utf8','strict')
    split = stripped.split()
    print(split)
    # get output and now convert into the necessary values to return to pygmo
    stripped2 = output2.stdout.strip().decode('utf8','strict')
    split2 = stripped2.split()
    print(split2)

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
            fp1xdim = (float(split[i])*1000)
            fp1xdim_bool = False
        if fp2xdim_bool:
            fp2xdim = (float(split[i])*1000)
            fp2xdim_bool = False
        if fp3xdim_bool:
            fp3xdim = (float(split[i])*1000)
            fp3xdim_bool = False
        if beamspotsize_bool:
            beamspotsize = power(float(split[i]),0.5)
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
            fp2res = (float(split2[i])*1000)
            fp2res_bool = False
        if fp3res_bool:
            fp3res = (float(split2[i])*1000)
            fp3res_bool = False
        if fp1res_bool:
            fp1res = (float(split2[i])*1000)
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

    # scale factor to account for the beam spot issue
    #   even the nominal setting is outside the bounds...
    scale = 1e9 
    max_width = 0
    # setup value to be returned, here 4 different objectives
    objs = 5
    resol = zeros(objs) 
    print(qs)
    for i in range(len(magnet_dims)):
        # if no x-ydim values, just return outside constraints (all 1e9)
        if len(xdim) < len(magnet_dims) or len(ydim) < len(magnet_dims):
            resol = zeros(objs)+1e9         
            break            
        if fp1xdim == 0 or fp2xdim == 0 or fp3xdim == 0 or isnan(fp1xdim) or isnan(fp2xdim) or isnan(fp3xdim) or xdim2[i]==0 or ydim2[i]==0:
            resol = zeros(objs)+1e9         
            break            
        # find xbound, ybound
        xbound = (abs(xdim[i])*1000)
        ybound = (abs(ydim[i])*1000)
        # if ratio of *bound to magnet_dim is larger than max_width, set max_width
        max_width = max(xbound/magnet_dims[i][0],ybound/magnet_dims[i][1],max_width)
#        print(i, magnet_names[i], "{:.2f}".format(xbound), magnet_dims[i][0], "{:.2f}".format(ybound), magnet_dims[i][1])
        # if *bound more than 4x magnet_dim, or is nan, return outside constraints
        if xbound > magnet_dims[i][0] * scale or ybound > magnet_dims[i][1] * scale or isnan(xbound) or isnan(ybound):
#            print("{:.2f}".format(xbound), magnet_dims[i][0], "{:.2f}".format(ybound), magnet_dims[i][1])
            resol = zeros(objs)+1e9         
            break
        # if within constraints, set resol temporarily
        try:
            resol = [fp1xdim/fp1res,fp2xdim/fp2res,fp3xdim/fp3res,max_width,beamspotsize]
        except:
            resol = zeros(objs)+1e9         
    print(resol)
    # if within constraints, set resol as a ratio to nominal
    if max(resol)/min(fNom) < 1e9:
        for i in range(len(resol)):
            # make sure we are working with positive numbers
            if resol[i] > 0: 
                resol[i] = float(resol[i])
            else:
                resol = zeros(objs)+1e9         
                break
            # we want to minimize max_width and beamspotsize, so just take resol/fNom
            resol[i] = resol[i]/fNom[i]
    else:
        resol = zeros(objs)+1e9         
    print(resol)            

    # remove old cosy fox and lis file
    commands.run(['rm','-f',cosyFilename])
    commands.run(['rm','-f',lisFilename])
    commands.run(['rm','-f',cosyFilename2])
    commands.run(['rm','-f',lisFilename2])

    # return the objective values
    return (resol)

def find_point(h5, index_i):
    df = pd.read_hdf(h5)
    qNew = power(zeros(19)+2,df.iloc[index_i,:19])
    return qNew 

if __name__ == "__main__":
    # if running from console, just run the nominal setting
    print(cosyrun(qNom))
    # if running from console, just run the assigned setting
#    qNew = find_point('results_280/best280.h5',52)
#    print(cosyrun(qNew))
    qNew = qNom
    PROFILES_PATH = "./"
    plot_i = 1
    write_fox(qNew, str(plot_i), PROFILES_PATH , 'SECAR_pg_Optics_draw.fox')
#    write_fox(qNew, str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE_draw.fox')
#    write_fox(qNew, str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_pg_Optics_DE2.fox')


