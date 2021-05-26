#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python3.9

import sys, math
import os, shutil, signal
#import commands
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

seed = 56448189
# SGA hyperparameters
generations = 300
cr_p = 0.9 # probability of crossover, 0.9 by default
mu_p = 0.7 # probability of mutation, 0.02 by default
mu_str = "gaussian" # mutation strategy, polynomial by default
mu_param_m = 0.05
magnet_dim = 11 
neighbors = 20
outputFile = "output_4f_moead_FP2_FP3_{}_init.csv".format(generations)


def main(pop_init=None):
    # Removing old files
#    x, v = createBees()   # Creating bees
    startTime = timeit.default_timer()

#    alg = pg.algorithm(pg.nsga2(gen=generations, cr=cr_p, m=mu_p))
#    alg = pg.nsga2(gen=generations, cr=cr_p, m=mu_p, seed=seed)
    alg = pg.moead(gen=generations,neighbours=neighbors,seed=seed)
#    alg = pg.algorithm(pg.sga(gen=generations,cr=cr_p,m=mu_p,mutation=mu_str,param_m=mu_param_m))
#    b = pg.bfe()
#    alg.set_bfe(b)
#    alg.set_verbosity(1)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim, outputFile))
#    print(p_optimizeRes)
#    pop = pg.population(p_optimizeRes, size=5)
#    pop = alg.evolve(pop)
#    print(pop.champion_f)
#    print(pop)
    n_islands = 1 
#    top = pg.topology(pg.fully_connected(n_islands,1.0))
    pop_n = 84 
    pop_new = None
    if (pop_init==None):    
        archi = pg.archipelago(n=n_islands,algo=alg,prob=p_optimizeRes,pop_size=pop_n)#,t=top
        print("initialize pop")                                                    #      
    else:                                                                          #      
        archi = pg.archipelago(n=n_islands,algo=alg,pop=pop_init)                  #,t=top
        pop_new = pop_init
        print("provided pop")
    archi.evolve()
    archi.wait_check()
    
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

# here's the 'magic'
    popi = None
#    for i, island in enumerate(archi): # iterate through islands
#        a = island.get_algorithm()      # get algorithm from island
#        popi = island.get_population()      # get algorithm from island
#        data = {}
#        if i==0 and pop_init==None:
#            pop_new = popi
#            for j in range(len(popi.get_f())):
#                data = {"id": [j+i*len(popi.get_f())],"island_#": [i], "f0": popi.get_f()[j][0], "f1": popi.get_f()[j][1], "f2": popi.get_f()[j][2], "f3": popi.get_f()[j][3]}
#                for k in range(len(popi.get_x()[j])):
#                    data["x"+str(k)] = popi.get_x()[j][k]
#                df = pd.DataFrame(data)
#                tot_df = tot_df.append(df)
#                tot_df.to_csv("test{0}_{1}_{2}_{3}_{4}.csv".format(generations,cr_p,mu_p,mu_str,mu_param_m),index=False)
#            continue
#        print("here's the result population: ")
#        print(popi)
#        for j in range(len(popi.get_f())):
#            pop_new.push_back(popi.get_x()[j], popi.get_f()[j])
#            data = {"id": [j+i*len(popi.get_f())],"island_#": [i], "f0": popi.get_f()[j][0], "f1": popi.get_f()[j][1], "f2": popi.get_f()[j][2], "f3": popi.get_f()[j][3]}
#            for k in range(len(popi.get_x()[j])):
#                data["x"+str(k)] = popi.get_x()[j][k]
#            df = pd.DataFrame(data)
#            tot_df = tot_df.append(df)
#        tot_df.to_csv("test{0}_{1}_{2}_{3}.csv".format(generations,cr_p,mu_p,mu_param_m),index=False)
            
#        np.concatenate((arr,popi.get_f()),axis=0)
#        print(arr)
#        tot_df = pd.concat([tot_df,df], axis='index', ignore_index=True) # merge with total df
#       uda = a.extract(pg.sga)         # extract algorithm from algorithm object
#       log = uda.get_log()             # get the log. Comes as list of tuples
#    
#       print(log)
       # reshape log
#       df = pd.DataFrame(np.asarray(log), columns = ["Gen", "F.evals.", "Best fit","mutation", "crossing over", "Variant", "dx", "df"])
#       df["island_#"] = i              # add island ID
    
#    uda = alg.extract(pg.sea)
#    log = uda.get_log()
#    print(log)
#    plt.semilogy([entry[0] for entry in log],[np.abs(entry[2]) for entry in log], 'k--') 
#    plt.savefig("test") 

    pg.mp_island.shutdown_pool()
    pg.mp_bfe.shutdown_pool()

    return pop_new

def init_pop(dim,pop_size):
    qNom = np.zeros(dim)
    pop = pg.population(prob=optimizeRes(dim,out=outputFile))
    pop.push_back(x=qNom,f=[0.,0.,0.,0.])#[0:dim])
    print(pop)
    print(pop.problem.get_fevals())
    pop.set_x(0,qNom)
    print(pop)
    for i in range(pop_size-1):
        pop.push_back(pop.random_decision_vector())  
    return pop

def read_pop(filename):
    df = pd.read_csv(filename)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
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
#    print(pop)
    return pop    
if __name__=='__main__':

    popi = init_pop(magnet_dim,84)
#    popi = read_pop("running_pop.csv")
#    popi = read_pop("f2_out.csv")
#    popi = read_pop("init_pop.csv")
    pop2 = main()
#    print("here's pop2: ")
#    print(pop2)
#    pop3 = main(pop2)
#    print("here's pop3: ")
#    print(pop3)
