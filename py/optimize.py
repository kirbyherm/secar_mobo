#!/usr/bin/env python3

#!/mnt/simulations/secarml/soft/anaconda3/bin/python
#!/usr/bin/python3.7
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

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

from cosy import cosyrun 
import pygmo as pg
from problem import optimizeRes
import pandas as pd

#set important directories
WORK_DIR = '../'
OUTPUT_DIR = WORK_DIR + 'output/' 
SCRATCH_DIR = '/scratch/hermanse'

# take batch_no as input from command line (this specifies outputFile)
script, batch_no = sys.argv

# set a specific seed for the algorithm
#   note that each island uses a random init pop
#   thus none of the islands are identical
#   (unless init pop is identical)
seed = 56448180

# MOEAD hyperparameters
#   default parameters have worked well
generations = 1000
cr_p = 1.0 # crossover parameter, 1.0 by default
f_p = 0.5 # diff evolution operator parameter, 0.5 by default
eta_m = 20 # distribution index used by the polynomial mutation, 20 by default
realb = 0.9 # chance that the neighborhood is considered at each generation, rather than the whole population, 0.9 by default
neighbors = 20 # size of the weight's neighborhood, 20 by default
limit = 2 # max number of copies reinserted in the population if preserve_diversity=True, 2 by default
preserve_diversity=True # activates diversity preservation mechanisms

# specify number of magnets to tune
magnet_dim = 19 
# specify output file
fileName = "output_5f_moead_gen_{}_{}.csv".format(generations, batch_no)
outputFile = "{}/{}".format(SCRATCH_DIR, fileName)

# main function
def main(pop_init=None):

    # start timer
    startTime = timeit.default_timer()

    # initialize algorithm with hyperparameters
    alg = pg.moead_gen(gen=generations,neighbours=neighbors,seed=seed)
#    alg = pg.nsga2(gen=generations)

    # there should be a way to run batch_fitness evaluations, haven't gotten it to work on NSCL
    b = pg.bfe()
    NUMBER_OF_PROCESSES=10
#    b.resize_pool(NUMBER_OF_PROCESSES)
    alg.set_bfe(b)
    #alg.set_verbosity(1)
    
    # initialize problem
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim, outputFile))

    # relic from old evolutions, which launched all islands internally
    #   this can allow for interconnectivity between islands
    #   now each run of this script calls its own island
    n_islands = 1 
    # if more than one island and want to exchange info need to define a topology to connect them    
    #top = pg.topology(pg.fully_connected(n_islands,1.0))

    # when running 5 objectives, pop needed to be 70    
    pop_n = 1001 #84 
    if p_optimizeRes.get_nobj() == 4:
        # with 4 objs need pop=84
        pop_n = 84 

    # check if we are using an input population or need to initialize one
    pop_new = None
    if (pop_init==None):    
        # randomly create a new population of size pop_n
#        pop_new = pg.population(p_optimizeRes, size=pop_n) 
        pop_new = pg.population(p_optimizeRes, size=pop_n, b=b) 
        print("initialize pop")
    else:
        pop_new = pop_init
        print("provided pop")

    # create archipelago from algorithm (of a single island)
    archi = pg.algorithm(alg)
    # evolve archipelago
    archi.evolve(pop_new)

    # check total time    
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

    # when I tried to use multiprocessing I needed this
    #pg.mp_island.shutdown_pool()
    #pg.mp_bfe.shutdown_pool()

    return pop_new

# if want to initialize a pop with a single nominal gene and random others
#   the evolution will quickly "forget" the nominal though as it prefers to find new solutions
def init_pop(dim,pop_size):
    qNom = np.zeros(dim)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    pop = pg.population(prob=optimizeRes(dim,out=outputFile))
    pop.push_back(x=qNom,f=np.zeros(p_optimizeRes.get_nobj()))#[0:dim])
    print(pop)
    print(pop.problem.get_fevals())
    pop.set_x(0,qNom)
    print(pop)
    for i in range(pop_size-1):
        pop.push_back(pop.random_decision_vector())  
    return pop


# if want to read in a specific population (as from a previous run)
#   read in from csv file, with header x0-x10, f0-f3
#     for a 11 magnet, 4 objective case
def read_pop(filename):
    df = pd.read_csv(filename)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim,out=outputFile))
    pop = pg.population(p_optimizeRes)
    nrow, ncol = df.shape
    print(df)
    for i in range(nrow):
        xs = []
        for j in range(magnet_dim):
            xs.append(df["x"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(ncol-magnet_dim):
            fs.append(df["f"+str(j)][i])
        pop.push_back(xs,f=fs)
    return pop    


# when called from console (as the batch script will)
if __name__=='__main__':


    isExist = os.path.exists(SCRATCH_DIR)
    if not isExist:
        os.makedirs(SCRATCH_DIR)

    # if using a preexisting or want to pass an identical population
    #popi = init_pop(magnet_dim,84)
    #popi = read_pop("init_pop.csv")

    # else just initialize and evolve random pop
    pop2 = main()

    shutil.move(outputFile, "{}{}".format(OUTPUT_DIR, fileName))
#    os.rmdir(SCRATCH_DIR)
