#!/usr/bin/env python3

import secar_utils as secar_utils
import make_profiles, draw_cluster_inverse
import draw_cluster
import draw_full, analyze_db, plot_tsne
from cosy import cosyrun

# Import libraries,
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import draw_cluster
import matplotlib
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
from itertools import repeat

from multiprocessing import Pool, Process, Queue, current_process, freeze_support

warnings.filterwarnings( "ignore", module = "matplotlib*" )

colormap = matplotlib.colormaps['PuOr']
plt.rcParams["figure.figsize"] = [14, 6]
colors = list(plt.colormaps['tab20'].colors)

configs = secar_utils.load_configs()

kclusters = configs['clusters']
n_obj = configs['n_obj']
n_con = configs['n_con']
magnet_dim = configs['magnet_dim']
objectives = configs['objectives']
fNom = configs['fNominal']

# set pandas view options to print everything
#pd.set_option("max_rows", None)
#pd.set_option("max_columns", None)

def make_db_row( quads, resol, cluster, columns ):
    new_row = {}
    new_row['kcluster'] = cluster
    new_row['closest'] = False
    new_row['ssobjs'] = 10
    for j in range(len(columns)):
        if j < magnet_dim:
            new_row[columns[j]] = quads[j]
        elif j < magnet_dim+n_obj+n_con:
            new_row[columns[j]] = resol[j-magnet_dim]
        else:
            break
    return new_row

def make_array_row( quads, resol, cluster, columns ):
    new_row = []
    for j in range(len(columns)):
        if j < magnet_dim:
            new_row.append(quads[j])
        elif j < magnet_dim+n_obj+n_con:
            new_row.append(resol[j-magnet_dim])
        else:
            break
    new_row.append(np.sum(np.power(new_row[magnet_dim:magnet_dim+n_obj+n_con],2)))
    new_row.append(cluster)
    new_row.append(0)
    return new_row

def parallel_eval(i,data, pca, mu, nComp, cluster):

    scale_check = False
    xrand = np.zeros(shape=(1,magnet_dim))
    while not scale_check:
        xrand = np.random.normal(scale=np.sqrt(pca.explained_variance_), size=(1, magnet_dim))
        Xhat = np.dot(xrand[:,:nComp], pca.components_[:nComp,:])
        Xhat += mu
        Xhat = Xhat.reshape((magnet_dim,))
        if min(Xhat) >= 0:
            scale_check = True
        if not scale_check:
            print("magnet factors outside bounds, re randomizing \n")        
    
    resol = np.zeros(configs['n_obj']+configs['n_con'])
    resol = cosyrun(Xhat, fNom)
    new_row = make_array_row( Xhat, resol, cluster, data.columns.insert(magnet_dim+n_obj,"FP4_BeamSpot") )
#    dataPCAtemp = pd.DataFrame(data=new_row,index=[0])
#    dataPCA = pd.concat([dataPCA, dataPCAtemp],ignore_index=True)
    return new_row

def main(results_no=0, samples=2, n_threads=1):


    magnets = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','q16','q17','q18','q19']
    magnet_scale_factors = []
    for i in range(len(magnets)):
#        if magnets[i] == 'q10':
#            magnet_scale_factors.append([0.125,8])
#        else: 
#            magnet_scale_factors.append([0.25,4])
        magnet_scale_factors.append([np.power(2,configs['lower_bounds'][i]),np.power(2,configs['upper_bounds'][i])])
    filepath = './results_{}/'.format(results_no)
    filename = filepath + 'best{}.h5'.format(results_no)
    data = pd.read_hdf(filename)
#    data.iloc[:,:19] = data.iloc[:,:19].apply(lambda x: np.power(2,x))   # Converting log factors to linear scale
    
    #data.kcluster.unique()  # number of clusters
    #data[data['kcluster'] == 1].kcluster.unique()

    pca = None
    mu = None
    nComp = 5
    
    dataPCA = pd.DataFrame(columns = data.columns.insert(magnet_dim+n_obj,"FP4_BeamSpot"))
    dataPCAT = pd.DataFrame(columns = data.columns.insert(magnet_dim+n_obj,"FP4_BeamSpot"))
    n_clusters = kclusters
#    fig, ax = plt.subplots(2,1)
    plot_combos = [[0,1],[2,3],[3,4],[2,4],[4,6],[6,9],[9,14],[12,14]]
    (fig, subplots) = plt.subplots(3, 3, figsize=(8, 8))
    for cluster in range(n_clusters):
        # Cluster 1:
        x1 = data[data.kcluster == cluster].iloc[:,:19]
        #x1 = data.iloc[:,:19]
        
        # Compute mean
        X = StandardScaler(with_std = False) 
        X = X.fit(x1)
        mu = X.mean_
    #    print(X.mean_, X.var_)
        X = X.transform(x1)
        
        # Compute PCA
        pca = PCA(n_components = magnet_dim).fit(X)
        # set up and run parallel minimization for n_sims samples
#        dataPCA = pd.read_hdf('results_2011/best2011PCA.h5')
        NUMBER_OF_PROCESSES = n_threads
        pool = Pool(processes = NUMBER_OF_PROCESSES)
        r = pool.starmap_async(parallel_eval, zip(range(samples), repeat(data), repeat(pca), repeat(mu), repeat(nComp), repeat(cluster)))
        #r = pool.apply_async(parallel_eval, args=( range(samples), data, cluster))
        r.wait()
#        r_return = np.array(r.get())
        for i in r.get():
            print(i)
            cols = dataPCA.columns
            new_row = {}
            for j in range(len(i)):
                new_row[dataPCA.columns[j]]=i[j]
            dataPCAtemp = pd.DataFrame(data=new_row,index=[0])
            dataPCA = pd.concat([dataPCA, dataPCAtemp],ignore_index=True)
        plot_x, plot_y = 0, 0
        pca = PCA(n_components = magnet_dim).fit(X)
        for pair in plot_combos:
        
            ax = subplots[plot_y][plot_x]
            ax.plot(dataPCA.loc[dataPCA["kcluster"]==cluster].iloc[:,pair[0]], dataPCA.loc[dataPCA["kcluster"]==cluster].iloc[:,pair[1]], 'o',color=colors[cluster])#,markersize=1.0)
            for comp in range(nComp):
                comps = [pca.components_[comp,pair[0]],pca.components_[comp,pair[1]]]
                comps = comps/np.linalg.norm(comps)*pca.explained_variance_[comp]
                if 9 in pair:
                    comps *= 4
                ax.arrow(mu[pair[0]],mu[pair[1]],comps[0],comps[1],zorder=10)
#            ax.plot(dataPCAT.loc[dataPCAT["kcluster"]==cluster].iloc[:,pair[0]], dataPCAT.loc[dataPCAT["kcluster"]==cluster].iloc[:,pair[1]], 'o',color=colors[cluster])#,markersize=1.0)
            plot_x += 1
            if plot_x > 2:
                plot_y+=1
                plot_x=0    
#        ax.set_xlim(0,2.0)
#        ax.set_ylim(0,2.0)
#            ax.plot(data.loc[data["kcluster"]==cluster].iloc[:,pair[0]], data.loc[data["kcluster"]==cluster].iloc[:,pair[1]], 'o',markersize=1.0,color=colors[cluster])#, fillstyle='none')
            ax.set_xlim(magnet_scale_factors[pair[0]][0],magnet_scale_factors[pair[0]][1])
            ax.set_ylim(magnet_scale_factors[pair[1]][0],magnet_scale_factors[pair[1]][1])
            ax.set_title("{1}:{0}".format(magnets[pair[0]], magnets[pair[1]]))
#        ax.set_xlim(0,2.0)
#        ax.set_ylim(0,2.0)
    plt.savefig(filepath+'pca_transformed')
    #    cosyrun(x1.iloc[0,:])
    db_out = 'results_{}/best{}PCA.h5'.format(results_no,results_no)
    dataPCA.to_hdf(db_out,key='df')
    query_txt= ''
    max_obj=configs['max_obj']
    for i in range(len(objectives)):
        query_txt += objectives[i] + "<{}".format(max_obj)
        if i < len(objectives)-1:
            query_txt+="&"
    query_txt+="&FP4_BeamSpot==1.0"
    print(dataPCA)
    dataPCA = dataPCA.query(query_txt)
    dataPCA = dataPCA.drop("FP4_BeamSpot",axis=1)
    print(dataPCA)
#    draw_cluster.main(data, dataPCA, 'results_{}/best{}PCA.h5'.format(results_no, results_no)) 
    draw_cluster_inverse.main(db_out, db_out) 
    dataPCA.to_hdf('../output/secar_{}d_db_{}s_PCA.h5'.format(n_obj, results_no),key='df')
#    print(x1.iloc[0,:],X[0,:],Xhat[0,:])
    draw_pca(pca, x1, results_no)
    return     

def draw_pca(pca, x1, results_no):
    
    # Singular values and cumulative energy
    S = pca.singular_values_
    sumS = [sum(S[:r]) for r in range(len(S)+1)]/sum(S[:])
    
    fig, ax = plt.subplots(1, 2)
    ax[0].plot(S)
    ax[0].set_yscale("log")
    ax[0].set_title("Singular values")
    
    ax[1].plot(sumS)
    ax[1].set_title("Cumulative values")
    
    ticks = np.arange(0, 20, 1)
    ax[1].set_xticks(ticks)
    plt.grid(True, which = 'both')
#    plt.show()
    plt.savefig('results_{}/best{}PCA.h5_Comps.png'.format(results_no, results_no) )

    plt.cla()
    fig, ax = plt.subplots(1, 1)
    # Plot PCA components # Remember to include the mean if using these components!!!
    for i in range(5):
        plt.plot(np.arange(1,20,1), pca.components_[i], label = str(i+1), linewidth=5-5*sumS[i])
    plt.legend()
    ticks = np.arange(0, len(pca.components_), 1)
    plt.xticks(ticks)
    plt.grid(True, which = 'both')
    plt.ylabel('Coefficient')
    plt.xlabel('Magnet')
    plt.savefig('results_{}/best{}PCA.h5_Mags.png'.format(results_no, results_no) )
#    plt.show()
    
    # Plot correlation matrix
    #print(x1.corr() )
    fig, ax = plt.subplots()
    
    threshold = 0.5
    x_corr = x1.corr()
    x_corr = x_corr[abs(x_corr) > threshold].apply(pd.to_numeric, errors='coerce').fillna(0, downcast='infer')
    plt.pcolormesh(x_corr,cmap=colormap)
    
    ax.set_xticks(np.arange(len(pca.components_))+0.5)
    ax.set_xticklabels(np.arange(len(pca.components_))+1)
    ax.set_yticks(np.arange(len(pca.components_))+0.5)
    ax.set_yticklabels(np.arange(len(pca.components_))+1)
    ax.set_ylabel('q')
    ax.set_xlabel('q')
    
    ax.set_title('correlation matrix for q1-q19')
    plt.colorbar()
    plt.savefig('results_{}/best{}PCA.h5_Corrs.png'.format(results_no, results_no) )
#    plt.show()
    return

if __name__=='__main__':

    import sys
    script, first, second, third = sys.argv
    results_no = int(first) 
    samples = int(second)
    n_threads = int(third)
    main(results_no, samples, n_threads)
