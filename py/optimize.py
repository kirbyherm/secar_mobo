#!/usr/bin/env python3

#import commands
import sys, math
import os, shutil, signal
import subprocess as commands
import re
import random
import numpy as np
import matplotlib.pyplot as plt
import time
import itertools
import timeit
import pygmo as pg
import pandas as pd

from cosy import cosyrun 
from problem import optimizeRes
#from problem_pca import optimizeRes
import secar_utils as secar_utils

configs = secar_utils.load_configs()
fox_name = configs['fox_name']
n_obj = configs['n_obj']
n_con = configs['n_con']
objectives_constraints = configs['objectives'] + configs['constraints']

# take batch_no as input from command line (this specifies the output db)
script, batch_no, use_prev = sys.argv

#set important directories
WORK_DIR = '../'
OUTPUT_DIR = WORK_DIR + 'output/' 
SCRATCH_DIR = configs['scratch_dir'] + batch_no+"/"
FOX_DIR = WORK_DIR + 'fox/' 

on_fireside = secar_utils.check_fireside()

TEMP_DIR = OUTPUT_DIR
TEMP_FOX_DIR = FOX_DIR

# set a specific seed for the algorithm
#   note that each island uses a random init pop
#   thus none of the islands are identical
#   (unless init pop is identical)
seed = 56448180

# MOEAD hyperparameters
#   default parameters have worked well
generations = 200
cr_p = 1.0 # crossover parameter, 1.0 by default
f_p = 0.5 # diff evolution operator parameter, 0.5 by default
eta_m = 1 # distribution index used by the polynomial mutation, 20 by default, lower value = higher variation, see sigma in eq. 7 of the moead/nsga2 paper
realb = 0.9 # chance that the neighborhood is considered at each generation, rather than the whole population, 0.9 by default
neighbors = 20 # size of the weight's neighborhood, 20 by default
limit = 2 # max number of copies reinserted in the population if preserve_diversity=True, 2 by default
preserve_diversity=True # activates diversity preservation mechanisms
weight_generation='low discrepancy'
pop_n = 1000

# specify number of magnets to tune
magnet_dim = configs['magnet_dim']
quads = []
for i in range(magnet_dim):
    quads.append("q{}".format(i+1))
columns = quads
for i in range(len(objectives_constraints)):
    columns.append(objectives_constraints[i])

# write population to h5 file
def save_pop(pop, db_out):

    xs = pop.get_x()
    fs = pop.get_f()
    df = pd.DataFrame(np.hstack((xs, fs)), columns=columns)            

    # Check whether the specified path exists or not
    isExist = os.path.exists(db_out)
    df_orig = None    
    if isExist:
        df_orig = pd.read_hdf(db_out)      
        df = pd.concat([df_orig, df], axis=0)
    df.to_hdf(db_out,key='df')
    
    return    

# main function
def main(db_out, pop_init=None ):

    # start timer
    startTime = timeit.default_timer()

    # initialize problem
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim, TEMP_FOX_DIR))

#    pop_n = 1001

    save_freq = int(generations *0.10)
    save_freq = int(generations)
    perc_prog = 0
    i = 0
    pop_new = None
#    while i < int(generations):
#        if perc_prog > 50:
#            save_freq = max(int(save_freq/2),2)
#        if perc_prog > 70:
#            save_freq = max(int(save_freq/2),2)
    # initialize algorithm with hyperparameters
    alg_constructor = None
    if 'secar_moead_custom' in sys.executable:
        alg_constructor = pg.moead_gen(gen=save_freq,neighbours=neighbors,weight_generation=weight_generation,eta_m=eta_m,seed=seed,outfile=SCRATCH_DIR+"output{}.csv".format(batch_no))
    else:
        alg_constructor = pg.moead_gen(gen=save_freq,neighbours=neighbors,weight_generation=weight_generation,seed=seed)
    b = pg.bfe()
    alg_constructor.set_bfe(b)
    alg = pg.algorithm(alg_constructor)
    #alg.set_verbosity(1)

    if i == 0:
        # check if we are using an input population or need to initialize one
        if (pop_init==None):    
            # randomly create a new population of size pop_n
            pop_new = pg.population(p_optimizeRes, size=pop_n, b=b) 
            print("initialize pop")
        else:
            pop_new = pop_init
            print("provided pop")

    pop_new = alg.evolve(pop_new)
    save_pop(pop_new, db_out)
#        i += save_freq
#        if int(i / generations *100 ) > perc_prog:
#            print ( i * save_freq )
#            perc_prog = int(i/generations * 100)
    # check total time    
#    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))
    
#    uda = alg.extract(pg.moead_gen)
#    print(uda.get_log())
    return pop_new

# deprecated function.....
# if want to initialize a pop with a single nominal gene and random others
#   the evolution will quickly "forget" the nominal though as it prefers to find new solutions
def init_pop(dim,pop_size):
    qNom = np.zeros(dim)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    pop = pg.population(prob=optimizeRes(dim))
    pop.push_back(x=qNom,f=np.zeros(p_optimizeRes.get_nobj()))
    print(pop)
    print(pop.problem.get_fevals())
    pop.set_x(0,qNom)
    print(pop)
    for i in range(pop_size-1):
        pop.push_back(pop.random_decision_vector())  
    return pop


# if want to read in a specific population (as from a previous run)
#   read in from h5 file, with specified objectives
def read_pop(filename, fox_dir=TEMP_FOX_DIR,pop_n=100):
    df = pd.read_hdf(filename)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim, fox_dir))
    pop = pg.population(p_optimizeRes)
    nrow, ncol = df.shape
    print(df)
    for i in range(nrow-pop_n,nrow):
        xs = []
        for j in range(magnet_dim):
            xs.append(df["q"+str(j+1)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(ncol-magnet_dim):
            fs.append(df[objectives_constraints[j]][i])
        pop.push_back(xs,f=fs)
    return pop    

# when called from console (as the batch script will)
if __name__=='__main__':


    on_fireside = secar_utils.check_fireside()

    db_out = "secar_{}d_db_{}s.h5".format(n_obj, batch_no)


    if on_fireside:
        isExist = os.path.exists(SCRATCH_DIR)
        if isExist:
            shutil.rmtree(SCRATCH_DIR, ignore_errors=True)
        os.makedirs(SCRATCH_DIR)
        TEMP_FOX_DIR = SCRATCH_DIR +'fox/'
        shutil.copytree(FOX_DIR, TEMP_FOX_DIR)
        TEMP_DIR = SCRATCH_DIR
    
    print(TEMP_FOX_DIR, SCRATCH_DIR)
#    TEMP_DIR = OUTPUT_DIR

#    if on_fireside:
#        isExist = os.path.exists(SCRATCH_DIR)
#        if isExist:
#            shutil.rmtree(SCRATCH_DIR)
#        os.makedirs(SCRATCH_DIR)
#        shutil.copytree(FOX_DIR, SCRATCH_DIR+'/fox/')
#        TEMP_DIR = SCRATCH_DIR

    # if using a preexisting or want to pass an identical population
    pop2 = None
    if int(use_prev) > 0:
        popi = read_pop(OUTPUT_DIR + "secar_{}d_db_{}s.h5".format(n_obj+n_con, int(batch_no)-1), TEMP_FOX_DIR, pop_n)
        pop2 = main(TEMP_DIR + db_out, popi)
    # else just initialize and evolve random pop
    else:
        pop2 = main(TEMP_DIR + db_out)
#    save_pop(pop2, db_out)

    if on_fireside:
        shutil.copyfile(SCRATCH_DIR+db_out, OUTPUT_DIR + db_out)
        if 'secar_moead_custom' in sys.executable:
            shutil.copyfile(SCRATCH_DIR+"output{}.csv".format(batch_no), OUTPUT_DIR + "output{}.csv".format(batch_no))
        shutil.rmtree(SCRATCH_DIR)


