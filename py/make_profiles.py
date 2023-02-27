#!/usr/bin/env python3

#import commands
import sys, math
import os, shutil, signal
import subprocess as commands
import re
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import time
import itertools
import timeit
import pygmo as pg


from cosy import cosyrun, write_fox
from problem import optimizeRes
from secar_utils import run_kmeans, is_pareto_efficient_simple, load_configs

configs = load_configs()
kclusters = configs['clusters']

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
plt.rcParams.update({
    "text.usetex": True,
})

optimized_params = configs['n_obj']
fNom = np.zeros(optimized_params)+1
fNames = [r"{FP1-res}${}^{-1}$",r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
magnet_dim = configs['magnet_dim']
objectives = configs['objectives']
fox_name = configs['fox_name']

fNom = np.array([2384.9360856494263, 109.61548781662407, 510.8029152516118, 1.6251646888022029, 0.12574090408565933])

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
    if 'PCA' in filename:
        costs = df[objectives]
        costs = np.array(costs)
        pareto = is_pareto_efficient_simple(costs)
        # add pareto column to df
        df['pareto'] = pareto
        # restrict df to only those points on the pareto front
        df = (df.loc[(df['pareto']==True)])
        df = df.reset_index(drop=True)
        df = df.drop(columns=['kcluster','closest'])
        df = run_kmeans(df, magnet_dim, 4) 
    return pop, df

# write cosy draw files from the best points
def output_4d_cosy(popi,filename,df):

    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    ref_point = np.zeros(optimized_params)+1e10 
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    magnet_dim = len(popi.get_x()[0])
    ndf_champ = []
    sorted_ndf = []
    sorted_pop = []
    sorted_xs = []
    sort_param = 4
    for i in ndf[0]:
        if np.all(np.array(popi.get_f()[i]) < 1) == True or True:
            ndf_champ.append(i)
    ndf = [ndf_champ]
    for i in ndf[0]:
        go_next = False
        if len(sorted_ndf) > 1:
            for j in range(len(sorted_ndf)):
                if np.array_equal(popi.get_f()[i],sorted_pop[j]) and np.array_equal(popi.get_x()[i],sorted_xs[j]):
                    go_next=True
                    break
        if go_next:
            continue            
        if i == ndf[0][0]:
            sorted_ndf.append(i)
            sorted_pop.append(popi.get_f()[i])
            sorted_xs.append(popi.get_x()[i])
        else:
            for j in range(len(sorted_ndf)):
                if j == len(sorted_ndf)-1:
                    sorted_ndf.append(i)
                    sorted_pop.append(popi.get_f()[i])
                    sorted_xs.append(popi.get_x()[i])
                    break
                elif j == 0 and popi.get_f()[i][sort_param] < popi.get_f()[sorted_ndf[j]][sort_param]:
                    sorted_ndf.insert(j,i)
                    sorted_pop.insert(j,popi.get_f()[i])
                    sorted_xs.insert(j,popi.get_x()[i])
                    break
                elif (popi.get_f()[i][sort_param] < popi.get_f()[sorted_ndf[j]][sort_param]) and j>0:
                    if(popi.get_f()[i][sort_param] >= popi.get_f()[sorted_ndf[j-1]][sort_param]):
#                        print(popi.get_f()[i][0],popi.get_f()[sorted_ndf[j]][0],popi.get_f()[sorted_ndf[j-1]][0])
                        sorted_ndf.insert(j,i)
                        sorted_pop.insert(j,popi.get_f()[i])
                        sorted_xs.insert(j,popi.get_x()[i])
                    break 
#    print(ndf[0], sorted_ndf)
    
    df_closest = df.loc[df['closest']==True]
#    print(df_closest, len(ndf[0]), len(sorted_ndf))
    results_path = os.path.split(filename)[0]
    PROFILES_PATH = results_path+'/profiles/'
    if os.path.exists(PROFILES_PATH) and os.path.isdir(PROFILES_PATH):
        files_in_dir = os.listdir(PROFILES_PATH)     # get list of files in the directory
        for file in files_in_dir:                  # loop to delete each file in folder
            os.remove(f'{PROFILES_PATH}/{file}')     
        shutil.rmtree(PROFILES_PATH)
    os.mkdir(PROFILES_PATH)
    qNew = np.zeros(19)+1
    write_fox((qNew), 0, PROFILES_PATH, '{}_draw.fox'.format(fox_name) )
    write_fox((qNew), str(0)+"_DE", PROFILES_PATH, '{}_DE_draw.fox'.format(fox_name) )
    count_dups = 0
    plot_i = 1
    for i in df_closest.index:
        i = i + 1
        j = sorted_ndf[i-1] 
#        print(popi.get_f()[j])
#        if i > 1:
#            for k in range(i-1):
#                if np.array_equal(sorted_pop[i-1],sorted_pop[k]):
#                    count_dups += 1
#                    break
        write_fox(popi.get_x()[j], plot_i, PROFILES_PATH, '{}_draw.fox'.format(fox_name))
        write_fox(popi.get_x()[j], str(plot_i)+"_DE", PROFILES_PATH, '{}_DE_draw.fox'.format(fox_name))
        plot_i += 1
    return

def main(filename):

    print("\nWriting cosy fox files for the {} identified cluster centroids\n".format(kclusters))
    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
    df = None
    if file_extension == ".h5":
        popi, df = read_pop_df(filename)
        output_4d_cosy(popi, filename, df)

if __name__=='__main__':
    script, filename = sys.argv
    main(filename)



