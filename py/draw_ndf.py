#!/usr/bin/env python3

#!/mnt/simulations/secarml/soft/anaconda3/bin/python
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

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

from cosy_draw import cosyrun, write_fox
import pygmo as pg
from problem import optimizeRes
import pandas as pd

os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'

plt.rcParams.update({
    "text.usetex": True,
})

# SGA hyperparameters
generations = 30 
cr_p = 0.9 # probability of crossover, 0.9 by default
mu_p = 0.7 # probability of mutation, 0.02 by default

script, filename = sys.argv
optimized_params = 4
fNom = np.zeros(optimized_params)+1
#fNames = [r"{FP2-res}${}^{-1}$",r"FP2-e-xangle",r"{FP3-res}${}^{-1}$",r"FP3-e-xangle",r"MaxBeamWidth"]
fNames = [r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
print(len(fNames))

def read_pop(filename,pop=None):
#    df = pd.read_csv(filename,names=["x0","x1","x2","x3","x4","x5","x6","f0","f1","f2","f3"])
    magnet_dim = len(pd.read_csv(filename).columns)-optimized_params
    p_optimizeRes = pg.problem(optimizeRes(magnet_dim))
    obj_dim = p_optimizeRes.get_nobj()
    quads = []
    for i in range(magnet_dim):
        quads.append("x{}".format(i))
    columns = quads
    for i in range(obj_dim):
        columns.append("f{}".format(i))
    df = pd.read_csv(filename,names=columns)
#    df["f0"] = df['f0'] * 1.0 / 4.419411469795324
#    df['f1'] = df['f1'] * 1.0 / 2.7701491204695756
#    print(magnet_dim)
    if pop == None:
        pop = pg.population(p_optimizeRes)
#    df_2 = df.loc[(df["f0"] < 10) & (df["f1"] < 10) & (df["f2"] < 10) & (df["f3"] < 10)]
#    df = df_2
    nrow, ncol = df.shape
#    print(df)
    for i in df.index:
#        if i in [0,102,76]:
#            continue
#    for i in range(nrow):
        xs = []
        for j in range(magnet_dim):
            xs.append(df["x"+str(j)][i]) 
        xs = np.asarray(xs)
        fs = []
#        print(ncol, magnet_dim)
        for j in range(ncol-magnet_dim):
#            if j > 0:
#                fs.append(df["f"+str(j)][i]/fNom[j])
#            else:
            if df["f"+str(j)][i] < 0 or np.isnan(df["f"+str(j)][i]):
                print(df["f"+str(j)][i], i)
                df["f"+str(j)][i] = 1e10
            fs.append(df["f"+str(j)][i])
        if np.all(np.array(fs) < 1) == True or True:
            pop.push_back(xs,f=fs)
#    print(pop)
    return pop    

def plot_2d_evo(popi):

    magnet_dim = len(popi.get_x()[0])
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
        ax.set_ylabel('xangle_e_min')
        ax.set_ylim(1e-3,1000)
        ax.set_yscale('log')
#    print(ndf_champ, ndf[0])
#    print(ndf)
        plt.savefig("popi{}".format(j))
    return

def output_2d_cosy(popi,filename):

    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    magnet_dim = len(popi.get_x()[0])
    ndf_champ = []
    sorted_ndf = []
    for i in ndf[0]:
        if i == ndf[0][0]:
            sorted_ndf.append(i)
        else:
            for j in range(len(sorted_ndf)):
                if j == len(sorted_ndf)-1:
                    sorted_ndf.append(i)
                    break
                elif j == 0 and popi.get_f()[i][0] < popi.get_f()[sorted_ndf[j]][0]:
                    sorted_ndf.insert(j,i)
                    break
                elif (popi.get_f()[i][0] < popi.get_f()[sorted_ndf[j]][0]) and j>0:
                    if(popi.get_f()[i][0] >= popi.get_f()[sorted_ndf[j-1]][0]):
#                        print(popi.get_f()[i][0],popi.get_f()[sorted_ndf[j]][0],popi.get_f()[sorted_ndf[j-1]][0])
                        sorted_ndf.insert(j,i)
                    break 
    print(ndf[0], sorted_ndf)
    
    for i in range(len(sorted_ndf)):
        j = sorted_ndf[i] 
        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), i, "2f_FP3/")
    return
def output_4d_cosy(popi,filename):

    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    magnet_dim = len(popi.get_x()[0])
    ndf_champ = []
    sorted_ndf = []
    for i in ndf[0]:
        if np.all(np.array(popi.get_f()[i]) < 1) == True:
            print(popi.get_f()[i])
            ndf_champ.append(i)
    ndf = [ndf_champ]
    for i in ndf[0]:
        if i == ndf[0][0]:
            sorted_ndf.append(i)
        else:
            for j in range(len(sorted_ndf)):
                if j == len(sorted_ndf)-1:
                    sorted_ndf.append(i)
                    break
                elif j == 0 and popi.get_f()[i][0] < popi.get_f()[sorted_ndf[j]][0]:
                    sorted_ndf.insert(j,i)
                    break
                elif (popi.get_f()[i][0] < popi.get_f()[sorted_ndf[j]][0]) and j>0:
                    if(popi.get_f()[i][0] >= popi.get_f()[sorted_ndf[j-1]][0]):
#                        print(popi.get_f()[i][0],popi.get_f()[sorted_ndf[j]][0],popi.get_f()[sorted_ndf[j-1]][0])
                        sorted_ndf.insert(j,i)
                    break 
    print(ndf[0], sorted_ndf)
    
    write_fox(np.power(np.zeros(magnet_dim)+2,np.zeros(magnet_dim)), 0, "4f_FP2_FP3/")
    for i in range(1,len(sorted_ndf)+1):
        j = sorted_ndf[i-1] 
        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), i, "4f_FP2_FP3/")
    return


def plot_2d(popi,filename):

    magnet_dim = len(popi.get_x()[0])
    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
#    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []
#    for j in range(2,200):
#    print(pg.select_best_N_mo(popi.get_f(),10))
    best_10 = (pg.select_best_N_mo(popi.get_f(),10))
    for i in best_10:
        print(i, np.power(np.zeros(magnet_dim)+2,popi.get_x()[i]), popi.get_f()[i])
    x = np.linspace(pg.ideal(popi.get_f())[0],0)
    y = np.zeros(50)+(1)
    plot_x, plot_y = 0,0
    objs = len(popi.get_f()[0])
    fig, axs = plt.subplots(objs-1,sharex=True)
    fig.suptitle('ratios to Nom vs Resolution')
    log_res, log_x = [], []
    for i in ndf[0]:
        log_res.append(np.log(popi.get_f()[i][0])/np.log(10))
        log_x.append(np.log(popi.get_f()[i][1])/np.log(10))
        if max(popi.get_f()[i]) < 1000000:
            print(i, np.power(np.zeros(magnet_dim)+2,popi.get_x()[i]), popi.get_f()[i])
    if objs == 2:
        j = 1
        axs.set_xlabel('resolution')
        axs.axvline(x=fNom[0],linestyle="dashed",color="red")
        axs.axvline(x=best_point[0],linestyle="dotted",color="blue")
        axs.axhline(y=1.0,linestyle="dashed",color="red")
        ndf_champ.append([popi.get_f()[i] for i in ndf[0]])
        pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,j],axes=axs)
#    ax.plot(color="C{}".format(j))
#        ax = pg.plot_non_dominated_fronts(popi.get_f()[0:j])

        axs.set_ylabel(fNames[2])
#        axs.set_ylim(0.9,1.1)
        axs.set_yscale('log')
#        axs.set_xlim(0.8,1.2)
        axs.set_xscale('log')
    else:
        axs[objs-2].set_xlabel('resolution')
        for j in range(1,2):
    #        axs[plot_y].plot(x,y,linestyle="dashed",color="red")
            axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
            axs[plot_y].axvline(x=best_point[0],linestyle="dotted",color="blue")
            axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
            ndf_champ.append([popi.get_f()[i] for i in ndf[0]])
            pg.plot_non_dominated_fronts(ndf_champ[0],comp=[0,j],axes=axs[plot_y])
    #    ax.plot(color="C{}".format(j))
    #        ax = pg.plot_non_dominated_fronts(popi.get_f()[0:j])
    
            axs[plot_y].set_ylabel(fNames[2])
    #        axs[plot_y].set_ylim(1e-1,1e1)
            axs[plot_y].set_yscale('log')
    #        axs[plot_y].set_xlim(0.1,10.0)
    #        axs[plot_y].set_xscale('log')
            plot_y += 1
#    print(ndf_champ, ndf[0])
#    print(ndf)
    fig.tight_layout()
    fig.savefig(filename+"_ndf.png")
    plt.clf()
#    plt.plot(log_res, log_x,"o", linestyle="none", )
#    fig.savefig(filename+"_logndf.png")
    return

def plot_4d(popi,filename):

    good_results=0
    magnet_dim = len(popi.get_x()[0])
    hv = pg.hypervolume(popi)
#    ref_point = hv.refpoint()
    ref_point = np.zeros(optimized_params)+1e10 
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
#    print(pg.sort_population_mo(popi.get_f()))
    ndf_champ = []
#    for j in range(2,200):
#    print(pg.select_best_N_mo(popi.get_f(),10))
    best_10 = (pg.select_best_N_mo(popi.get_f(),10))
#    for i in best_10:
#        print(i, np.power(np.zeros(magnet_dim)+2,popi.get_x()[i]), popi.get_f()[i])
    x = np.linspace(pg.ideal(popi.get_f())[0],0)
    y = np.zeros(50)+(1)
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at FP4')
    axs[optimized_params-2].set_xlabel(fNames[3])
#    for i in ndf[0]:
#        if max(popi.get_f()[i]) < 10:
#            print(i, np.power(np.zeros(magnet_dim)+2,popi.get_x()[i]), popi.get_f()[i])
    reduced_ndf = []
    first = True
    for j in range(0,optimized_params-1):
#        axs[plot_y].plot(x,y,linestyle="dashed",color="red")
        axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
#        axs[plot_y].axvline(x=best_point[0],linestyle="dotted",color="blue")
        axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
        for i in ndf[0]:
            check_val=1e9
            if np.all(np.array(popi.get_f()[i]) < check_val) == True:
#            if np.all(np.array(popi.get_f()[i]) < 1) == True or first:
                good_results+=1
#                print(filename[-6:-4], good_results, popi.get_f()[i])
                ndf_champ.append(popi.get_f()[i])
                reduced_ndf.append(i)
#                first = False
        try:
            pg.plot_non_dominated_fronts(ndf_champ,comp=[3,j],axes=axs[plot_y])
            if j == 1:
                print(filename[-6:-4], len(ndf_champ))
        except:
            print(filename[-6:-4], "no better than nominal solutions")
            return
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
    plt.cla()
    fig2, axs2 = plt.subplots(3,4)
    plot_x, plot_y = 0,0
    reduced_qs = np.array(popi.get_x())[reduced_ndf]
#    print(reduced_qs)
    df = pd.DataFrame(reduced_qs, columns = ['y0','y1','y2','y3','y4','y5','y6','y7','y8','y9','y10'])
    qNom = np.zeros(magnet_dim)
    for i in range(magnet_dim):
        if plot_x > 3:
            plot_x, plot_y = 0, plot_y+1 
        axs2[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs2[plot_y,plot_x],bins=20,range=(-1,1))
        axs2[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
#        axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
        axs2[plot_y,plot_x].axes.yaxis.set_visible(False)
        axs2[plot_y,plot_x].axes.set_xlim(-1,1)
        axs2[plot_y,plot_x].set_title("q{0}".format(i+1))
        plot_x += 1
    
    fig2.delaxes(axs2[plot_y,plot_x])
    fig2.tight_layout()
    plt.savefig(filename + ".png")
    return

def main(filename):
    popi = read_pop(filename)
#    filename = '/mnt/simulations/secarml/secar_mobo/output/output_4f_moead_FP2_FP3_150_'
#    for i in range(40,50):
#        popi = read_pop(filename + '{}.csv'.format(i),popi)
#    print(popi)
#    popi = read_pop("init_pop.csv")
    hv = pg.hypervolume(popi)
#    ref_point = hv.refpoint()
    ref_point = np.zeros(optimized_params)+1e10 
#    print(ref_point)
    hv.compute(ref_point = ref_point)
#    print(hv.compute(ref_point = ref_point))
#    print(hv.greatest_contributor(ref_point))
#    print(popi.get_f()[hv.greatest_contributor(ref_point)])
#    plt.show()
#    plot_2d(popi,filename)
#    output_2d_cosy(popi, filename)
    plot_4d(popi,filename)
#    output_4d_cosy(popi, filename)
#    plot_2d(popi)

if __name__=='__main__':
    main(filename)

