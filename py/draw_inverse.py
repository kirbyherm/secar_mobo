#!/usr/bin/env python3

#!/mnt/simulations/secarml/soft/anaconda3/bin/python
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
import matplotlib
import time
import itertools
import timeit
from cosy import cosyrun, write_fox
import pygmo as pg
from problem import optimizeRes
import pandas as pd
from kmeans_cluster import run_kmeans
from view_db import is_pareto_efficient_simple

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
plt.rcParams.update({
    "text.usetex": True,
})

script, filename = sys.argv
optimized_params = 5
fNom = np.zeros(optimized_params)+1
fNames = [r"{FP1-res}${}^{-1}$",r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
#print(len(fNames))
magnet_dim = 19
objectives = ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','DSSD_BeamSpot']
fNom1 = np.array([2384.9360856494263, 109.61548781662407, 510.8029152516118, 1.6251646888022029, 0.12574090408565933])

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

# 2d functions are not updated, this one was for creating pngs for a movie of the evolution
#   movie made from pngs using ffmpeg
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

# outdated
def output_2d_cosy(popi,filename):

    hv = pg.hypervolume(popi)
    ref_point = hv.refpoint()
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    magnet_dim = len(popi.get_x()[0])
    ndf_champ = []
    sorted_ndf = []
    sort_param = 3
    for i in ndf[0]:
        if i == ndf[0][0]:
            sorted_ndf.append(i)
        else:
            for j in range(len(sorted_ndf)):
                if j == len(sorted_ndf)-1:
                    sorted_ndf.append(i)
                    break
                elif j == 0 and popi.get_f()[i][sort_param] < popi.get_f()[sorted_ndf[j]][sort_param]:
                    sorted_ndf.insert(j,i)
                    break
                elif (popi.get_f()[i][sort_param] < popi.get_f()[sorted_ndf[j]][sort_param]) and j>0:
                    if(popi.get_f()[i][sort_param] >= popi.get_f()[sorted_ndf[j-1]][sort_param]):
#                        print(popi.get_f()[i][0],popi.get_f()[sorted_ndf[j]][0],popi.get_f()[sorted_ndf[j-1]][0])
                        sorted_ndf.insert(j,i)
                    break 
    print(ndf[0], sorted_ndf)
    
    for i in range(len(sorted_ndf)):
        j = sorted_ndf[i] 
        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), i, "2f_FP3/")
    return

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
#            print(popi.get_f()[i])
            ndf_champ.append(i)
    ndf = [ndf_champ]
    for i in ndf[0]:
#        print(i)
        go_next = False
        if len(sorted_ndf) > 1:
            for j in range(len(sorted_ndf)):
#                print( j, sorted_ndf)
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
    print(df_closest, len(ndf[0]), len(sorted_ndf))
    results_path = os.path.split(filename)[0]
    PROFILES_PATH = results_path+'/profiles/'
    if os.path.exists(PROFILES_PATH) and os.path.isdir(PROFILES_PATH):
        files_in_dir = os.listdir(PROFILES_PATH)     # get list of files in the directory
        for file in files_in_dir:                  # loop to delete each file in folder
            os.remove(f'{PROFILES_PATH}/{file}')     
        shutil.rmtree(PROFILES_PATH)
    os.mkdir(PROFILES_PATH)
    qNew = np.array([1.0371301857113335,1.4897519431921593,0.5402003843384104,0.6080163749223835,0.5965351874518491,0.5279178522813484,0.8474952322221544,0.8290931192132953,0.7350223112146984,0.5049139345530922,0.969681779928563 ,0.8465270119223961,0.7261232553654523,0.6805787940919176,0.6772214286022437,1.6737045402403927,1.3151418622198896,0.8914897696929639,0.6144362243855045])
    qNew = np.zeros(19)+1
#    qNew = np.array([0.82281,0.83116,0.94706,0.75348,0.64679,1.23591,2.51027,0.80262,0.8017,2.70642,1.20897,0.62487,0.38632,1.11944,0.86569,1.19812,1.64921,1.13609,0.28077])
    write_fox((qNew), 0, PROFILES_PATH, 'SECAR_an_Optics_draw.fox' )
    write_fox((qNew), str(0)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE_draw.fox' )
#    write_fox((qNew), str(0)+"_DE_FP1", PROFILES_PATH, 'SEC_neutrons_WF_off_DE_rays_v1_draw_FP1.fox' )
    count_dups = 0
#    for i in range(1,len(sorted_ndf)+1):
#    print(df_closest,df_closest.index)
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
#        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), plot_i, PROFILES_PATH, 'SECAR_an_Optics_draw.fox')
#        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE_draw.fox')
        write_fox(popi.get_x()[j], plot_i, PROFILES_PATH, 'SECAR_an_Optics_draw.fox')
        write_fox(popi.get_x()[j], str(plot_i)+"_DE", PROFILES_PATH, 'SECAR_an_Optics_DE_draw.fox')
#        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), str(plot_i)+"_DE_FP1", PROFILES_PATH, 'SEC_neutrons_WF_off_DE_rays_v1_draw_FP1.fox')
        plot_i += 1
#    print(len(sorted_ndf), count_dups)
    return

def plot_4d(popi,filename,df):

    sort_param = 4
    good_results=0
    magnet_dim = len(popi.get_x()[0])
#    hv = pg.hypervolume(popi)
#    ref_point = np.zeros(optimized_params)+1e10 
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at FP4')
    axs[optimized_params-2].set_xlabel(fNames[sort_param])
    reduced_ndf = []
    first = True
    df_closest = df
    cmap = matplotlib.colors.ListedColormap(matplotlib.cm.get_cmap("Pastel1").colors[:3])

    for j in range(0,optimized_params-1):
        ndf_champ = []
        axs[plot_y].axvline(x=fNom1[sort_param],linestyle="dashed",color="red")
        axs[plot_y].axhline(y=fNom1[j],linestyle="dashed",color="red")
        for i in ndf[0]:
            check_val=1e9
            if np.all(np.array(popi.get_f()[i]) < check_val) == True:
                good_results+=1
#                print(filename[-6:-4], good_results, popi.get_f()[i])
                transform_f = popi.get_f()[i]
                for k in range(len(transform_f)):
                    if k < 3:
                        transform_f[k] = fNom1[k] / transform_f[k]
                    else:
                        transform_f[k] = transform_f[k] * fNom1[k]
                ndf_champ.append(transform_f)
                reduced_ndf.append(i)
        try:
            pg.plot_non_dominated_fronts(ndf_champ,comp=[sort_param,j],axes=axs[plot_y])
#            if j == 1:
#                print(filename[-6:-4], len(ndf_champ))
        except:
            print(filename[-6:-4], "no better than nominal solutions")
            return
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
            df_closest = df_closest.reset_index(drop=True)
#            print(df_closest.iloc[:,15:19])
#            for i_closest in df_closest.index:
#                axs[plot_y].text(df_closest.iloc[:,magnet_dim+sort_param][i_closest],df_closest.iloc[:,magnet_dim+j][i_closest],str(i_closest+1),color='red')
        axs[plot_y].set_ylabel(fNames[j])
#        axs[plot_y].set_yscale('log')
#        print(axs[plot_y].get_xlim(),axs[plot_y].get_ylim())
#        max_y = max(abs(np.log10(axs[plot_y].get_ylim())))
        ndf_champ_array = np.array(ndf_champ)
        max_lim = max(abs(ndf_champ_array[:,j].min()-fNom1[j]),abs(ndf_champ_array[:,j].max()-fNom1[j]))
        axs[plot_y].set_ylim(ndf_champ_array[:,j].min()-max_lim,ndf_champ_array[:,j].max()+max_lim)
#        axs[plot_y].set_xscale('log')
#        max_x = max(abs(np.log10(axs[plot_y].get_xlim())))
        max_lim = max(abs(ndf_champ_array[:,sort_param].min()-fNom1[sort_param]),abs(ndf_champ_array[:,sort_param].max()-fNom1[sort_param]))
        axs[plot_y].set_xlim(ndf_champ_array[:,sort_param].min()-max_lim,ndf_champ_array[:,sort_param].max()+max_lim)
#        axs[plot_y].set_xlim(np.power(10,-max_x),np.power(10,max_x))
#        print(math.ceil(max(axs[plot_y].get_ylim())),math.ceil(max(axs[plot_y].get_ylim())))
#        colors = np.zeros((int(math.ceil(max(axs[plot_y].get_ylim()))),int(math.ceil(max(axs[plot_y].get_xlim()))))) + 1
#        for colori in range(1,colors.shape[0]):
#            for colorj in range(1, colors.shape[1]):
#                colors[colori,colorj] = 0
#        colors[0,0] = 2
##        print(colors)
#        axs[plot_y].pcolormesh(colors, cmap=cmap)
        plot_y += 1

    fig.tight_layout()
    fig.savefig(filename+"_paretos_inverse.png")
    plt.cla()
    fig2, axs2 = plt.subplots(5,4)
    plot_x, plot_y = 0,0
    reduced_qs = np.array(popi.get_x())[reduced_ndf]
    ycolumns = []
    for i in range(magnet_dim):
        ycolumns.append('y{}'.format(i))
    df = pd.DataFrame(reduced_qs, columns = ycolumns)
    qNom = np.zeros(magnet_dim)
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
    for i in range(magnet_dim):
        if plot_x > 3:
            plot_x, plot_y = 0, plot_y+1 
        axs2[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs2[plot_y,plot_x],bins=100,range=(-3,3))
#        axs2[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
#        axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
        axs2[plot_y,plot_x].axes.yaxis.set_visible(False)
#        axs2[plot_y,plot_x].axes.set_yscale("log")
#        axs2[plot_y,plot_x].axes.set_ylim(0.1,20)
        xlower, xupper = popi.problem.get_bounds()
        xlower, xupper = np.min(xlower), np.max(xupper)
        axs2[plot_y,plot_x].axes.set_xlim(xlower,xupper)
        axs2[plot_y,plot_x].set_title("{0}".format(write_qnames[i]))
        y_min, y_max = axs2[plot_y,plot_x].get_ylim()
        df_closest['yplot'] = pd.Series(df_closest.index).apply(lambda x: x/len(df_closest.index)*(y_max-y_min)+y_min)
#        print(df_closest.iloc[:,:15])
#        if "closest" in df.columns:
        for i_closest in df_closest.index:
            axs2[plot_y,plot_x].text(df_closest["q{0}".format(i+1)][i_closest],df_closest['yplot'][i_closest],str(i_closest+1),color='red')
        
        plot_x += 1
    
    axs2[plot_y,plot_x].axis('off')
    axs2[plot_y,plot_x].text(0.3, 0.5, 'x-axis is in log2', horizontalalignment='center', verticalalignment='center', transform=axs2[plot_y,plot_x].transAxes)
    fig2.tight_layout()
#    plt.savefig(filename + "_magnet_hists.png")
    return

def plot_hists(df, df_reduce, filename):

    write_qnames = {'q1':'q1','q2':'q2','q3':'q3','q4':'q4','q5':'q5','q6':'q6','q7':'q7','q8':'q8','q9':'q9','q10':'q10','q11':'q11','q12':'q12','q13':'q13','q14':'q14','q15':'q15','q16':'h1','q17':'h2','q18':'h3','q19':'o1'}
    df = df.iloc[:,:19]
    df_reduce = df_reduce.iloc[:,:19]
    df = df.rename(columns=write_qnames)
#    fig = plt.figure(1)
    fig, axs = plt.subplots(5,4)
    print(axs)
    fig.delaxes(axs[4][3])
    axs = (fig.axes)
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
        print(axs[i].get_ylim())
        axs[i].set_ylim([100,50000])
        axs[i].set_ylabel('all points')
        if i != 9:
#            axs[i].set_xlim([-2,2])
            axs[i].set_xlim([0.25,4])
        else:
#            axs[i].set_xlim([-3,3])
            axs[i].set_xlim([0.125,8])
        axs2[i] = axs[i].twinx()
        axs2[i].set_ylim([1,1000])
        axs2[i].set_ylabel('optimized points',rotation=-90)
    hist = df_reduce.hist(bins=bins,log=True,ax=axs2,color='red',grid=False)
    for i in range(len(axs)):
        axs2[i].set_title("")

    print(hist)
    plt.savefig(filename+'_full_hist.png')
    
    return


def main(filename):

    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
    df = None
    if file_extension == ".h5":
        popi, df = read_pop_df(filename)
#        output_4d_cosy(popi, filename, df)
    else: 
        popi = read_pop(filename)
#    print(df)
    plot_4d(popi,filename,df)

if __name__=='__main__':
    main(filename)



