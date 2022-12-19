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

# specify Tex details for pretty plots
os.environ['PATH'] = os.environ['PATH'] + ':/mnt/misc/sw/indep/all/texlive/2013/bin/x86_64-linux/latex'
#os.environ['PATH'] = os.environ['PATH'] + ':/usr/bin/tex'
#plt.rcParams.update({
#    "text.usetex": True,
#})

# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

optimized_params = 5
fNom = np.array([1000/2384.9360856494263, 1000/109.61548781662407, 1000/510.8029152516118, 1.6251646888022029/1.6251646888022029, 0.12574090408565933 * 100])
fNames = [r"{FP1-res}${}^{-1}$",r"{FP2-res}${}^{-1}$",r"{FP3-res}${}^{-1}$",r"MaxBeamWidth",r"BeamSpotSize"]
fNames = fNames[:optimized_params]
#print(len(fNames))
magnet_dim = 19

def correct_obj_values(df, filename):

    
    df = df.loc[df['closest']==True]
    objs =  ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    nclusters = max(df['kcluster']+1) 
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
#            if j == 1:
#                print(filename[-6:-4], len(ndf_champ))
        except:
            print(filename[-6:-4], "no better than nominal solutions")
            return
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
            df_closest = df_closest.reset_index(drop=True)
#            print(df_closest.iloc[:,15:19])
            for i_closest in df_closest.index:
                axs[plot_y].text(df_closest.iloc[:,magnet_dim+sort_param][i_closest],df_closest.iloc[:,magnet_dim+j][i_closest],str(i_closest+1),color='red')
        axs[plot_y].set_ylabel(fNames[j])
        axs[plot_y].set_yscale('log')
#        print(axs[plot_y].get_xlim(),axs[plot_y].get_ylim())
        max_y = max(abs(np.log10(axs[plot_y].get_ylim())))
        axs[plot_y].set_ylim(np.power(10,-max_y),np.power(10,max_y))
        axs[plot_y].set_xscale('log')
        max_x = max(abs(np.log10(axs[plot_y].get_xlim())))
        axs[plot_y].set_xlim(np.power(10,-max_x),np.power(10,max_x))
#        print(math.ceil(max(axs[plot_y].get_ylim())),math.ceil(max(axs[plot_y].get_ylim())))
        colors = np.zeros((int(math.ceil(max(axs[plot_y].get_ylim()))),int(math.ceil(max(axs[plot_y].get_xlim()))))) + 1
        for colori in range(1,colors.shape[0]):
            for colorj in range(1, colors.shape[1]):
                colors[colori,colorj] = 0
        colors[0,0] = 2
#        print(colors)
        axs[plot_y].pcolormesh(colors, cmap=cmap)
        plot_y += 1

    fig.tight_layout()
#    fig.savefig(filename+"_paretos.png")
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

def main(df, df_compare, filename, PCA_run = False):
#    df_compare = correct_obj_values(df_compare,filename)
    sort_param = 4
    number_of_clusters = np.max(df['kcluster']+1)
    colors = np.array(list(plt.get_cmap('tab20').colors)).reshape(-1,3)
    print(colors, colors.shape, np.array(colors[0]).reshape(3))
    objectives = ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    j = 0
    plot_x, plot_y = 0,0
    fig, axs = plt.subplots(optimized_params-1,sharex=True)
    fig.suptitle('Pareto Fronts of each parameter vs. BeamSpotSize at DSSD')
    axs[3].set_xlabel("DSSD_BeamSpot")
    for obj in ['FP1_res','FP2_res','FP3_res','MaxBeamWidth']:
        df.plot(x='FP4_BeamSpot',y=obj,style='o',ax=axs[plot_y],markersize=3.0,legend=False)
        if "closest" in df.columns:
            df_closest = df.loc[df['closest']==True]
            df_closest = df_closest.reset_index(drop=True)
#            print(df_closest.iloc[:,15:19])
            df_closest.plot(x='FP4_BeamSpot',y=obj,style='o',ax=axs[plot_y],markersize=3.0,legend=False)
#            for i_closest in df_closest.index:
#                axs[plot_y].text(df_closest.iloc[:,magnet_dim+sort_param][i_closest],df_closest.iloc[:,magnet_dim+j][i_closest],str(i_closest+1),color='black')
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
#        plt.savefig(filename+'_'+obj+'_inverse.png')
        j += 1
        plot_y += 1

#    fig.tight_layout()
    axs[0].set_ylabel("FP1 Res.")
    axs[1].set_ylabel("FP2 Res.")
    axs[2].set_ylabel("FP3 Res.")
    axs[3].set_xlabel("DSSD Beamspot (cm)")
    plt.savefig(filename+'_inverse.png')


if __name__=='__main__':

    script, filename, filename_compare = sys.argv
    file_extension = os.path.splitext(filename)[-1]
    print(os.path.split(filename))
    popi = None
    print("reading {} file".format(file_extension))
#    if file_extension == ".h5":
#        popi, df_list = read_pop_df(filename)
#    for df in df_list:
#        plot_4d(popi,filename,df)

    df = pd.read_hdf(filename)
    df_compare = pd.read_hdf(filename_compare)
    PCA_run = False
    objectives = ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    for i in range(len(objectives)):
        if i < 3:
            df[objectives[i]] = df[objectives[i]].apply(lambda x: fNom[i] / x)
        else:
            if i != 3:
                df[objectives[i]] = df[objectives[i]].apply(lambda x: fNom[i] * x)
    print(df[objectives])
    
#    if 'PCA' in filename_compare:
#        query_txt = '' 
#        max_obj = 1
#        for i in range(len(objectives)):
#            query_txt += objectives[i] + "<{}".format(max_obj)
#            if i < len(objectives)-1:
#                query_txt+="&"
#        df_compare = df_compare.query(query_txt)
#        filename=filename_compare
#        PCA_run = True
    main(df, df_compare, filename, PCA_run)


