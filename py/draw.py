#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

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
import time
import itertools
import timeit
from cosy import cosyrun, write_fox
import pygmo as pg
from problem import optimizeRes
import pandas as pd

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
os.environ['PATH'] = os.environ['PATH'] + ':/opt/software/texlive/20210316/bin/x86_64-linux/tex'

plt.rcParams.update({
    "text.usetex": True,
})

script, filename = sys.argv
optimized_params = 4
fNom = np.zeros(optimized_params)+1
fNames = [r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
print(len(fNames))

# read pop from h5 file (i.e. after running view_db.py)
def read_pop_df(filename, pop=None):
    df = pd.read_hdf(filename)
    magnet_dim = 15
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
                print(df["f"+str(j)][i], i)
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
    ref_point = (1e10,1e10,1e10,1e10) 
    best_point = (popi.get_f()[hv.greatest_contributor(ref_point)])
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    magnet_dim = len(popi.get_x()[0])
    ndf_champ = []
    sorted_ndf = []
    sorted_pop = []
    sorted_xs = []
    sort_param = 3
    for i in ndf[0]:
        if np.all(np.array(popi.get_f()[i]) < 1) == True or True:
#            print(popi.get_f()[i])
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
    if not ( os.path.isdir("4f_FP2_FP3") ):
        os.mkdir("4f_FP2_FP3")
    write_fox(np.power(np.zeros(magnet_dim)+2,np.zeros(magnet_dim)), 0, "4f_FP2_FP3/", 'SEC_neutrons_WF_14m_v1_draw.fox' )
    count_dups = 0
#    for i in range(1,len(sorted_ndf)+1):
    print(df_closest,df_closest.index)
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
        write_fox(np.power(np.zeros(magnet_dim)+2,popi.get_x()[j]), plot_i, "4f_FP2_FP3/", 'SEC_neutrons_WF_14m_v1_draw.fox')
        plot_i += 1
#    print(len(sorted_ndf), count_dups)
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

def plot_4d(popi,filename,df):

    sort_param = 3
    good_results=0
    magnet_dim = len(popi.get_x()[0])
    hv = pg.hypervolume(popi)
    ref_point = np.zeros(optimized_params)+1e10 
    ndf, dl, dc, ndl = pg.fast_non_dominated_sorting(popi.get_f())
    ndf_champ = []
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at FP4')
    axs[optimized_params-2].set_xlabel(fNames[sort_param])
    reduced_ndf = []
    first = True
    df_closest = df
    for j in range(0,optimized_params-1):
        axs[plot_y].axvline(x=fNom[0],linestyle="dashed",color="red")
        axs[plot_y].axhline(y=1.0,linestyle="dashed",color="red")
        for i in ndf[0]:
            check_val=1e9
            if np.all(np.array(popi.get_f()[i]) < check_val) == True:
                good_results+=1
#                print(filename[-6:-4], good_results, popi.get_f()[i])
                ndf_champ.append(popi.get_f()[i])
                reduced_ndf.append(i)
        try:
            pg.plot_non_dominated_fronts(ndf_champ,comp=[sort_param,j],axes=axs[plot_y])
            if j == 1:
                print(filename[-6:-4], len(ndf_champ))
        except:
            print(filename[-6:-4], "no better than nominal solutions")
            return
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
#            print(df_closest.iloc[:,15:19])
            axs[plot_y].plot(df_closest.iloc[:,18],df_closest.iloc[:,15+j],'x',color='red')
        axs[plot_y].set_ylabel(fNames[j])
        axs[plot_y].set_yscale('log')
        axs[plot_y].set_xscale('log')
        plot_y += 1

    fig.tight_layout()
    fig.savefig(filename+"_ndf.png")
    plt.cla()
    fig2, axs2 = plt.subplots(4,4)
    plot_x, plot_y = 0,0
    reduced_qs = np.array(popi.get_x())[reduced_ndf]
    ycolumns = []
    for i in range(magnet_dim):
        ycolumns.append('y{}'.format(i))
    df = pd.DataFrame(reduced_qs, columns = ycolumns)
    qNom = np.zeros(magnet_dim)
    df_closest = df_closest.reset_index(drop=True)
    for i in range(magnet_dim):
        if plot_x > 3:
            plot_x, plot_y = 0, plot_y+1 
        axs2[plot_y,plot_x] = df['y{0}'.format(i)].plot.hist(ax=axs2[plot_y,plot_x],bins=100,range=(-1,1))
        axs2[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
#        axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
        axs2[plot_y,plot_x].axes.yaxis.set_visible(False)
        axs2[plot_y,plot_x].axes.set_xlim(-1,1)
        axs2[plot_y,plot_x].set_title("q{0}".format(i+1))
        y_min, y_max = axs2[plot_y,plot_x].get_ylim()
        df_closest['yplot'] = pd.Series(df_closest.index).apply(lambda x: x/15*(y_max-y_min)+y_min)
        print(df_closest.iloc[:,:15])
#        if "closest" in df.columns:
        axs2[plot_y,plot_x].plot(df_closest["q{0}".format(i+1)],df_closest['yplot'],'x',color='red')
        
        plot_x += 1
    
    fig2.delaxes(axs2[plot_y,plot_x])
    fig2.tight_layout()
    plt.savefig(filename + ".png")
    return

def main(filename):

    file_extension = os.path.splitext(filename)[-1]
    popi = None
    print(file_extension)
    df = None
    if file_extension == ".h5":
        popi, df = read_pop_df(filename)
        output_4d_cosy(popi, filename, df)
    else: 
        popi = read_pop(filename)
    print(df)
    plot_4d(popi,filename,df)

if __name__=='__main__':
    main(filename)



