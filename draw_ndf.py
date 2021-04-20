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
fNom = [-870.9084099103394, 0.01625617474106655, 0.02579740767554017, 0.0020116413429212]
fNom = [1,1,1,1]
fNames = ["resolution","xwidth_e","xangle_e","xangle_xwidth"]
qNom = np.array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )

script, filename = sys.argv

def read_pop(filename):
    df = pd.read_csv(filename,names=["x0","x1","x2","x3","x4","x5","x6","f0","f1","f2","f3"])
    p_optimizeRes = pg.problem(optimizeRes(7))
    pop = pg.population(p_optimizeRes)
#    df_2 = df.loc[(df["f0"] < 10) & (df["f1"] < 10) & (df["f2"] < 10) & (df["f3"] < 10)]
#    df = df_2
    nrow, ncol = df.shape
    print(df)
    for i in df.index:
#    for i in range(nrow):
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

    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []
#    for j in range(2,200):
    print(pg.select_best_N_mo(popi.get_f(),10))
    best_10 = (pg.select_best_N_mo(popi.get_f(),10))
    for i in best_10:
        print(i, np.power(np.zeros(7)+2,popi.get_x()[i]), popi.get_f()[i])
    x = np.linspace(pg.ideal(popi.get_f())[0],0)
    y = np.zeros(50)+(1)
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(3,sharex=True)
    fig.suptitle('ratios to Nom vs Resolution')
    axs[2].set_xlabel('resolution')
    for i in ndf[0]:
        if max(popi.get_f()[i]) < 10:
            print(i, np.power(np.zeros(7)+2,popi.get_x()[i]), popi.get_f()[i])
    for j in range(1,4):
#        axs[plot_y].plot(x,y,linestyle="dashed",color="red")
        axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
        axs[plot_y].axvline(x=best_point[0],linestyle="dotted",color="blue")
        axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
        ndf_champ.append([popi.get_f()[i] for i in ndf[0]])
        pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,j],axes=axs[plot_y])
#    ax.plot(color="C{}".format(j))
#        ax = pg.plot_non_dominated_fronts(popi.get_f()[0:j])

        axs[plot_y].set_ylabel(fNames[j])
#        axs[plot_y].set_ylim(1e-1,1e1)
        axs[plot_y].set_yscale('log')
#        axs[plot_y].set_xlim(0.1,10.0)
        axs[plot_y].set_xscale('log')
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
    print(popi.get_f()[hv.greatest_contributor(ref_point)])
#    plt.show()
    plot_4d(popi,filename)
#    plot_2d(popi)

if __name__=='__main__':
    main(filename)

