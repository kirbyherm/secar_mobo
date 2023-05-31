#!/usr/bin/env python3

import secar_utils as secar_utils

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
import pandas as pd
import pygmo as pg

from problem import optimizeRes

configs = secar_utils.load_configs()

kclusters = configs['clusters']
optimized_params = configs['n_obj']
objectives = configs['objectives']
fNom = configs['fNominal_plot']
magnet_dim = configs['magnet_dim']
max_obj = configs['max_obj']

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
#plt.rcParams.update({
#    "text.usetex": True,
#})

# set pandas view options to print everything
#pd.set_option("max_rows", None)
#pd.set_option("max_columns", None)

# deprecated function that was used to scale the obj values for different batches
def correct_obj_values(df, filename):

    
    df = df.loc[df['closest']==True]
    objs = objectives
    nclusters = kclusters
    x1 = np.array([0.4149030189991541, 1.024219663615905, 0.32276368444381337, 0.28454926622496207, 0.13051896098280003])
    x2 = np.append(x1,[0.40651582802980557, 0.9107875198036088, 0.24913824216441277, 3.9806676668358776, 3.808775729951597])
    x3 = np.append(x2,[0.4162521335795522, 0.9444542701708413, 0.5153166810864698, 0.19393722897375568, 0.1705774768981977])
    x4 = np.append(x3,[0.3945841068155405, 0.518947611092936, 0.5147924366850625, 0.22812565873404525, 0.22328122790935367])

    if '340' in filename:
        x1 = np.array([1.1602257895768526, 1.9832986250319096, 1.2019965432541717, 1.0225988358141749, 0.5612462612961021])
        x2 = np.append(x1,[4,4,4,4,4])
        x3 = np.append(x2,[1.1645946603141062, 1.7004600000167525, 1.7789270074156527, 0.6033383813418585, 0.6022615821623248])
        x4 = np.append(x3,[1.1390023263255333, 2.0992401393492774, 2.2332984461485816, 0.9315916411551571, 0.9715838648006931])

    for i in range(nclusters):
        for j in range(len(objs)):
            df.loc[df['kcluster']==i,objs[j]] = x4[i*len(objs)+j]
    return df
        

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)

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
            fs.append(df.iloc[i,j])
        if append:
            pop.push_back(xs,f=fs)
    return pop, df

def plot_pareto(df, df_compare, filename, PCA_run = False):
#    df_compare = correct_obj_values(df_compare,filename)
    sort_param = 4
    number_of_clusters = np.max(df['kcluster']+1)
    colors = np.array(list(plt.get_cmap('tab20').colors)).reshape(-1,3)
    print(colors, colors.shape, np.array(colors[0]).reshape(3))
    j = 0
    if optimized_params < 5:
        j+=1
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at DSSD')
    axs[optimized_params-2].set_xlabel("DSSD_BeamSpot")
    for obj in objectives[:-1]:
        df.plot(x='FP4_BeamSpot',y=obj,style='o',ax=axs[plot_y],markersize=3.0,legend=False)
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
            df_closest = df_closest.reset_index(drop=True)
            df_closest.plot(x=objectives[-1],y=obj,style='o',ax=axs[plot_y],markersize=3.0,legend=False)
        axs[plot_y].axes.set_ylabel(obj)
        axs[plot_y].axvline(x=fNom[-1],linestyle="dashed",color="red")
        axs[plot_y].axhline(y=fNom[j],linestyle="dashed",color="red")
        cmap = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("Pastel1").colors[:3])
        colors_x = np.zeros((2,2)) 
        print(axs[plot_y].get_xlim(), axs[plot_y].get_ylim())
        xlims = [axs[plot_y].get_xlim()[0], fNom[-1], axs[plot_y].get_xlim()[1]]
        ylims = [axs[plot_y].get_ylim()[0], fNom[j], axs[plot_y].get_ylim()[1]]
        for colori in range(colors_x.shape[0]):
            for colorj in range(colors_x.shape[1]):
                colors_x[colori,colorj] = 0
        if j < 3:
            colors_x[1,0] = 2
        else:
            colors_x[0,0] = 2
        print(colors_x)
        axs[plot_y].pcolormesh(xlims, ylims, colors_x, cmap=cmap)
        j += 1
        plot_y += 1

    axis_index = 0
    if optimized_params==5:
        axs[axis_index].set_ylabel("FP1 Res.")
        axis_index += 1
    else:
        axis_index=0
    axs[axis_index].set_ylabel("FP2 Res.")
    axis_index += 1
    axs[axis_index].set_ylabel("FP3 Res.")
    axis_index += 1
    axs[axis_index].set_xlabel("DSSD Beamspot (cm)")
    plt.savefig(filename+'_pareto.png')

def main(filename, filename_compare):

    print("\nDrawing the pareto front of the points\n")
    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
    df = pd.read_hdf(filename)
#    df_compare = pd.read_hdf(filename_compare)
    df_compare = pd.read_hdf(filename)
    PCA_run = False
    if 'PCA' in filename_compare:
        query_txt = '' 
        for i in range(len(objectives)):
            query_txt += objectives[i] + "<{}".format(max_obj)
            if i < len(objectives)-1:
                query_txt+="&"
        df_compare = df_compare.query(query_txt)
        df = df_compare
        filename=filename_compare
        PCA_run = True
        print(PCA_run)
    for i in range(len(objectives)):
        fNom_i = i
        if optimized_params < 5:
            fNom_i += 1
        if fNom_i < 3:
            df[objectives[i]] = df[objectives[i]].apply(lambda x: fNom[fNom_i] / x)
        else:
            df[objectives[i]] = df[objectives[i]].apply(lambda x: fNom[fNom_i] * x)
    
    plot_pareto(df, df_compare, filename, PCA_run)

if __name__=='__main__':

    script, filename, filename_compare = sys.argv
    main(filename, filename_compare)



