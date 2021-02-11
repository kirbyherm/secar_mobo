#!/usr/bin/env python

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


def main():
    # Removing old files
    
#    x, v = createBees()   # Creating bees
    startTime = timeit.default_timer()

#    alg = pg.algorithm(pg.nsga2())
    alg = pg.algorithm(pg.sga(gen=10))
    alg.set_verbosity(1)
    p_optimizeRes = pg.problem(optimizeRes(7))
#    print(p_optimizeRes)
#    pop = pg.population(p_optimizeRes, size=5)
#    pop = alg.evolve(pop)
#    print(pop.champion_f)
#    print(pop)
    archi = pg.archipelago(n=1,algo=alg,prob=p_optimizeRes,pop_size=10)
    archi.evolve()
    archi.wait_check()
    print(archi.get_migrants_db())
    print(archi.get_migration_log())
#    print(archi.get_champions_f())
#    print(archi.status)    
    
    print ('Running time (sec): %f' % (timeit.default_timer() - startTime))

#    uda = alg.extract(pg.sea)
#    log = uda.get_log()
#    print(log)
#    plt.semilogy([entry[0] for entry in log],[np.abs(entry[2]) for entry in log], 'k--') 
#    plt.savefig("test") 

    pg.mp_island.shutdown_pool()
    pg.mp_bfe.shutdown_pool()

    return 0

if __name__=='__main__':
    main()


