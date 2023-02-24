#!/usr/bin/env python3

#import commands
import sys, math
import os, shutil, signal
import subprocess as commands
import re
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import time
import itertools
import timeit
import pygmo as pg
import pandas as pd

from cosy import cosyrun, write_fox
from problem import optimizeRes
import secar_utils as secar_utils

configs = secar_utils.load_configs()
# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
plt.rcParams.update({
    "text.usetex": True,
})

kclusters = configs['clusters']
optimized_params = configs['n_obj']
objectives = configs['objectives']
fNom = configs['fNominal_plot']
magnet_dim = configs['magnet_dim']

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)
    magnet_dim = 19
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    nobj = p_optimizeRes.get_nobj()
    if pop == None:
        pop = pg.population(p_optimizeRes)
    nrow, ncol = df.shape
    for i in df.index:
        append=True
        xs = []
        for j in range(1,magnet_dim+1):
            xs.append(df["q"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(magnet_dim,magnet_dim+nobj):
#            if i == 0:
#                print(df.iloc[i,magnet_dim:magnet_dim+p_optimizeRes.get_nobj()])
            fs.append(df.iloc[i,j])
        if append:
            pop.push_back(xs,f=fs)
    return pop, df

# read in pop from csv output file
def read_pop(filename,pop=None):
    # get magnet dimensions 
    magnet_dim = len(pd.read_csv(filename).columns)-optimized_params
    # init problem for creating pop
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    obj_dim = p_optimizeRes.get_nobj()
    quads = []
    # construct columns to read from csv
    for i in range(magnet_dim):
        quads.append("x{}".format(i))
    columns = quads
    for i in range(obj_dim):
        columns.append("f{}".format(i))
    # read df from csv
    df = pd.read_csv(filename,names=columns)
    # initialize pop if none
    if pop == None:
        pop = pg.population(p_optimizeRes)

    # construct pop from df 
    nrow, ncol = df.shape
    for i in df.index:
        xs = []
        for j in range(magnet_dim):
            xs.append(df["x"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
        for j in range(ncol-magnet_dim):
            # if not a valid obj value throw out point
            if df["f"+str(j)][i] < 0 or np.isnan(df["f"+str(j)][i]):
                print("error: ", df["f"+str(j)][i], i)
                df["f"+str(j)][i] = 1e10
            fs.append(df["f"+str(j)][i])
        # add point to population
        pop.push_back(xs,f=fs)
    # return population to main()
    return pop    

def plot_hists(df, df_reduce, filename):

    write_qnames = {'q1':'q1','q2':'q2','q3':'q3','q4':'q4','q5':'q5','q6':'q6','q7':'q7','q8':'q8','q9':'q9','q10':'q10','q11':'q11','q12':'q12','q13':'q13','q14':'q14','q15':'q15','q16':'h1','q17':'h2','q18':'h3','q19':'o1'}
    for i in range(magnet_dim):
#        df.iloc[:,i] = df.iloc[:,i].apply(lambda x: np.log2(np.power(2,x)*scale_factor[i]))
        df.iloc[:,i] = df.iloc[:,i].apply(lambda x: np.power(2,x))
#    df_invalid = df.loc[np.sum(df,axis=1) >= 5e9]
    df_invalid = df
    print(df_invalid)
    df_invalid = df_invalid.iloc[:,:19]
    df = df.iloc[:,:19]
    df_clos = df_reduce.loc[df_reduce['closest']==True].sort_values(by=['FP4_BeamSpot']).reset_index(drop=True)
    df_best = df_reduce.copy()
    df_reduce = df_reduce.iloc[:,:19]
    df = df.rename(columns=write_qnames)
#    fig = plt.figure(1)
    fig, axs = plt.subplots(5,4)
#    print(axs)
    fig.delaxes(axs[4][3])
    axs = np.array(fig.axes)
    bins = np.linspace(-3,3,100)
    bins = np.logspace(-3,3,100,base=2.0)
    hists = df.hist(bins=bins,log=True,ax=axs,grid=False)
#    print(hists[0].bins)
#    fig.yscale('log')
    fig.tight_layout()
    fig.set_figheight(10)
    fig.set_figwidth(19)
#    ax = plt.gca()
#    print(ax.get_ylim())
    axs2 = axs
    for i in range(len(axs)):
#    for j in range(len(axs[i])):
#        print(axs[i].get_ylim())
        axs[i].set_ylim([2000,150000])
        axs[i].set_ylabel('all points')
        axs[i].set_xscale('log')
        if i != 9:
#            axs[i].set_xlim([-2,2])
            axs[i].set_xlim([0.25,4])
        else:
#            axs[i].set_xlim([-3,3])
            axs[i].set_xlim([0.125,8])
        axs2[i] = axs[i].twinx()
        axs2[i].set_ylim([2000,150000])
#        axs2[i].set_ylim([1,int(1e5)])
        axs2[i].set_ylabel('invalid points',rotation=-90)
    number_of_clusters = np.max(df_best['kcluster']+1)
    colors = list(plt.get_cmap('tab20').colors)
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','q16','q17','q18','q19']
#    print(df_best,axs2)
    df_clos['yplot'] = 0
#    for i in range(number_of_clusters):
##                ax = df.loc[(df['kcluster']==i)].plot(x='FP4_BeamSpot',y=obj,style='o',color=colors[i],label=df_clos.loc[df_clos['kcluster']==i].index[0]+1,markersize=3.0)
#        hist = df_best.loc[df_best['kcluster']==i][write_qnames].hist(bins=bins,log=True,ax=axs2,color=np.array([colors[i+1]]),label=df_clos.loc[(df_clos['kcluster']==i)].index[0]+1,grid=False)
#        for j in range(len(axs)):
#            y_min, y_max = axs2[j].get_ylim()
##            df_clos['yplot'][i] = df_clos.index*(y_max-y_min)+y_min
##            print((df_clos['kcluster'][i])*(np.logspace(1.0,3.0,num=10)[i]),(df_clos['kcluster'][i]),(np.logspace(0.0,2.5,num=10)[i]))
#            axs2[j].text(df_clos["q{0}".format(j+1)][i],(np.logspace(0.0,2.8,num=10)[df_clos['kcluster'][i]]),str(np.max(df_clos.loc[(df_clos['kcluster']==i)]['kcluster'])+1),color='black')
    hist = df_invalid[write_qnames].hist(bins=bins,log=True,ax=axs2,color='k',label="invalid",grid=False)
#    for j in range(len(axs)):
#        y_min, y_max = axs2[j].get_ylim()
#            df_clos['yplot'][i] = df_clos.index*(y_max-y_min)+y_min
#            print((df_clos['kcluster'][i])*(np.logspace(1.0,3.0,num=10)[i]),(df_clos['kcluster'][i]),(np.logspace(0.0,2.5,num=10)[i]))
#        axs2[j].text(df_invalid["q{0}".format(j+1)][i],(np.logspace(0.0,2.8,num=10)[df_clos['kcluster'][i]]),str(np.max(df_clos.loc[(df_clos['kcluster']==i)]['kcluster'])+1),color='black')
    for i in range(len(axs)):
        axs2[i].set_title("")

#    print(hist)
    plt.savefig(filename+'full_hist_invalid.png')
#    plt.savefig(filename+'full_hist.png')
    
    return


def main(filename,batch):

    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
    df = None
    if file_extension == ".h5":
        popi, df = read_pop_df(filename)
    df_full = pd.read_hdf('../output/secar_{}d_db_{}s.h5'.format(configs['n_obj'],batch))
    plot_hists(df_full, df, "results_{}/".format(batch))

if __name__=='__main__':
    script, filename, batch = sys.argv
    main(filename,batch)



