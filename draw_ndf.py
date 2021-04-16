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
import pandas as pd

# SGA hyperparameters
generations = 30 
cr_p = 0.9 # probability of crossover, 0.9 by default
mu_p = 0.7 # probability of mutation, 0.02 by default
magnet_dim = 7
fNom = [-445.6197096190245, 0.01622835298051074, 0.0244989752733823, 0.01640748422423699]
fNames = ["resolution","xwidth_e","xangle_e","xangle_xwidth"]

script, filename = sys.argv

def read_pop(filename):
    df = pd.read_csv(filename,names=["x0","x1","x2","x3","x4","x5","x6","f0","f1","f2","f3"])
    p_optimizeRes = pg.problem(optimizeRes(7))
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
            if j > 0:
                fs.append(df["f"+str(j)][i]/fNom[j])
            else:
                fs.append(df["f"+str(j)][i])
        pop.push_back(xs,f=fs)
#    print(pop)
    return pop    

def plot_2d(popi):

    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []
#    for j in range(2,200):
    print(pg.ideal(popi.get_f())[0])
    x = np.linspace(pg.ideal(popi.get_f())[0],0)
    y = np.zeros(50)+(1)
    for j in range(2,100):
        plt.cla()
        plt.plot(x,y,linestyle="dashed",color="red")
        print(ndf[0])
        ndf_champ.append([popi.get_f()[i] for i in ndf[j]])
#        ax = pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,j])
#    ax.plot(color="C{}".format(j))
        ax = pg.plot_non_dominated_fronts(popi.get_f()[0:j])
        ax.set_xlabel('resolution')
        ax.set_ylabel('xwidth_e_min')
        ax.set_ylim(1e-3,1000)
        ax.set_yscale('log')
#    print(ndf_champ, ndf[0])
#    print(ndf)
        plt.savefig("popi{}".format(j))
    return

def plot_4d(popi,filename):

    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []
#    for j in range(2,200):
    print(pg.ideal(popi.get_f())[0])
    x = np.linspace(pg.ideal(popi.get_f())[0],0)
    y = np.zeros(50)+(1)
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(3,sharex=True)
    fig.suptitle('ratios to Nom vs Resolution')
    axs[2].set_xlabel('resolution')
    for j in range(1,4):
#        axs[plot_y].plot(x,y,linestyle="dashed",color="red")
        axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
        axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
        print(ndf[0],plot_y)
        ndf_champ.append([popi.get_f()[i] for i in ndf[0]])
        pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,j],axes=axs[plot_y])
#    ax.plot(color="C{}".format(j))
#        ax = pg.plot_non_dominated_fronts(popi.get_f()[0:j])

        axs[plot_y].set_ylabel(fNames[j])
        axs[plot_y].set_ylim(1e-4,10)
        axs[plot_y].set_yscale('log')
        plot_y += 1
#    print(ndf_champ, ndf[0])
#    print(ndf)
    fig.tight_layout()
    fig.savefig(filename+"_ndf.png")
    return

def main(filename):
    popi = read_pop(filename)
#    popi = read_pop("init_pop.csv")
    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    print(ref_point)
#    hv.compute(ref_point = ref_point)
#    print(hv.compute(ref_point = ref_point))
    print(hv.greatest_contributor(ref_point))
    print(pg.ideal(popi.get_f()))
#    plt.show()
    plot_4d(popi,filename)
#    plot_2d(popi)

if __name__=='__main__':
    main(filename)

