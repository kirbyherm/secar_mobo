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
from problem import optimizeRes
import pandas as pd

os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'

plt.rcParams.update({
    "text.usetex": True,
})

script, filename = sys.argv
optimized_params = 4
# SGA hyperparameters
generations = 30 
cr_p = 0.9 # probability of crossover, 0.9 by default
mu_p = 0.7 # probability of mutation, 0.02 by default
fNom = [-870.9084099103394, 0.01625617474106655, 0.02579740767554017, 0.0020116413429212]
fNom = [1,1,1,1]
fNames = [r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]


def is_pareto_efficient_simple(costs):
    """
    Find the pareto-efficient points
    :param costs: An (n_points, n_costs) array
    :return: A (n_points, ) boolean array, indicating whether each point is Pareto efficient
    """
    is_efficient = np.ones(costs.shape[0], dtype = bool)
    for i, c in enumerate(costs):
#        print(i,c,is_efficient[i])
        if is_efficient[i]:
            is_efficient[is_efficient] = np.any(costs[is_efficient]<c, axis=1)  # Keep any point with a lower cost
            is_efficient[i] = True  # And keep self
    return is_efficient

def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)
    max_val = 1e9
    df = df.loc[(df['FP2_res'] < max_val) & (df['FP3_res'] < max_val) & (df['MaxBeamWidth'] < max_val) & (df['FP4_BeamSpot'] <max_val)]
    df = df.reindex()
    magnet_dim = len(df.columns)-optimized_params
#    df = df.loc[(df['FP2_res'] < 1.0) & (df['FP2_e_xangle'] < 1.0) & (df['FP3_res'] < 1.0) & (df['FP3_e_xangle'] <1.0)]
#    costs = df[['FP2_res','FP2_e_xangle','FP3_res','FP3_e_xangle']]
#    costs = np.array(costs)
#    pareto = is_pareto_efficient_simple(costs)
#    df['pareto'] = pareto
#    print(np.count_nonzero(pareto) )
#    df = (df.loc[(df['pareto']==True)])
#    df = df.sort_values(by='FP2_res',ignore_index=True)
#    df["f0"] = df['f0'] * 1.0 / 4.419411469795324
#    df['f1'] = df['f1'] * 1.0 / 2.7701491204695756
#    print(magnet_dim)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    if pop == None:
        pop = pg.population(p_optimizeRes)
#    df_2 = df.loc[(df["f0"] < 10) & (df["f1"] < 10) & (df["f2"] < 10) & (df["f3"] < 10)]
#    df = df_2
    nrow, ncol = df.shape
#    print(df.index)
    for i in df.index:
        if i >26000:
            break
#        print(i)
        append=True
        xs = []
        for j in range(1,magnet_dim+1):
            xs.append(df["q"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(magnet_dim,ncol+0):
#            if j > 0:
#                fs.append(df["f"+str(j)][i]/fNom[j])
#            else:
#            if df["f"+str(j)][i] >= 1.0:
#                append=False 
#            if df["f"+str(j)][i] < 0 or np.isnan(df["f"+str(j)][i]):
#                print(df["f"+str(j)][i], i)
            fs.append(df.iloc[i,j])
        if append:
            pop.push_back(xs,f=fs)
#    print(pop)
    return pop    

def read_pop(filename, pop=None):
    df = pd.read_csv(filename,names=["x0","x1","x2","x3","x4","x5","x6","f0","f1","f2","f3"])
    magnet_dim = len(df.columns)-optimized_params
    quads = []
    for i in range(magnet_dim):
        quads.append("x{}".format(i))
    columns = quads
    columns.append("f0")
    columns.append("f1")
    columns.append("f2")
    columns.append("f3")
    
#    df["f0"] = df['f0'] * 1.0 / 4.419411469795324
#    df['f1'] = df['f1'] * 1.0 / 2.7701491204695756
#    print(magnet_dim)
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    if pop == None:
        pop = pg.population(p_optimizeRes)
#    df_2 = df.loc[(df["f0"] < 10) & (df["f1"] < 10) & (df["f2"] < 10) & (df["f3"] < 10)]
#    df = df_2
    nrow, ncol = df.shape
#    print(df)
    for i in df.index:
        append=True
        xs = []
        for j in range(magnet_dim):
            xs.append(df["x"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(ncol-magnet_dim):
#            if j > 0:
#                fs.append(df["f"+str(j)][i]/fNom[j])
#            else:
#            if df["f"+str(j)][i] >= 1.0:
#                append=False 
            if df["f"+str(j)][i] < 0 or np.isnan(df["f"+str(j)][i]):
                print(df["f"+str(j)][i], i)
            fs.append(df["f"+str(j)][i])
        if append:
            pop.push_back(xs,f=fs)
#    print(pop)
    return pop    

def plot_2d_evo(popi):

    magnet_dim = len(popi.get_x()[0])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f()[:,:2])
#    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []

#    for j in range(2,200):
#    print(pg.ideal(popi.get_f())[0])
    print(ndf)
#    x = np.linspace(pg.ideal(popi.get_f())[0],0)
#    y = np.zeros(50)+(1)
    old_ndf, new_ndf = np.zeros(1), np.zeros(1) 
    n_plots = 0
    for j in range(2,len(popi.get_f())):
#        print(x,y)
#        plt.plot(x,y,linestyle="dashed",color="red")
#        print(ndf)
#        ndf_champ.append([popi.get_f()[i] for i in ndf[j]])
#        ax = pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,1])
#    ax.plot(color="C{}".format(j))
        ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f()[0:j,:2])
#        print(ndf)
        new_ndf = np.array(ndf[0])
        if not np.array_equal(new_ndf, old_ndf):
            if len(ndf[0]) <=1:
                continue
            n_plots += 1
            plt.cla()
            ax = pg.plot_non_dominated_fronts(popi.get_f()[ndf[0]],comp=[0,1])
            ax.set_ylabel(fNames[1])
            ax.set_xlabel(fNames[0])
            ax.set_ylim(1e-1,1e1)
            ax.set_yscale('log')
            ax.set_xlim(0.1,10.0)
            ax.set_xscale('log')
            ax.axvline(x=fNom[0],linestyle="dashed",color="red")
#            axs[plot_y].axvline(x=best_point[0],linestyle="dotted",color="blue")
            ax.axhline(y=1.0,linestyle="dashed",color="red")
#            ax.set_xlabel('resolution')
#            ax.set_ylabel('xangle_e_min')
#            ax.set_ylim(1e-3,1000)
#            ax.set_yscale('log')
#    pri    nt(ndf_champ, ndf[0])
#    pri    nt(ndf)
            plt.savefig("popi{}".format(n_plots))
        old_ndf = np.array(ndf[0])
    return

def main(filename):
    popi = read_pop_df(filename)
    hv = pg.hypervolume(popi)
    ref_point = (1e10,1e10,1e10,1e10) 
    hv.compute(ref_point = ref_point)
    print(hv.compute(ref_point = ref_point))
    plot_2d_evo(popi)

if __name__=='__main__':
    main(filename)

