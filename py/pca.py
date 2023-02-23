#!/usr/bin/env python3

# Import libraries,
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import draw_cluster
# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)


import matplotlib
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import warnings
warnings.filterwarnings( "ignore", module = "matplotlib*" )

from cosy import cosyrun

colormap = matplotlib.cm.get_cmap('PuOr')
plt.rcParams["figure.figsize"] = [14, 6]
colors = list(plt.get_cmap('tab20').colors)

results_no = 430

def make_db_row( quads, resol, cluster, columns ):
    new_row = {}
    new_row['kcluster'] = cluster
    new_row['closest'] = True
    new_row['ssobjs'] = 10
    for j in range(len(columns)):
        if j < 19:
            new_row[columns[j]] = quads[j]
        elif j < 19+5:
            new_row[columns[j]] = resol[j-19]
        else:
            break
    return new_row

def main():


    magnets = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','q16','q17','q18','q19']
    magnet_scale_factors = []
    for i in range(len(magnets)):
        if magnets[i] == 'q10':
            magnet_scale_factors.append([0.125,8])
        else: 
            magnet_scale_factors.append([0.25,4])
    filepath = './results_{}/'.format(results_no)
    filename = filepath + 'best{}.h5'.format(results_no)
    data = pd.read_hdf(filename)
#    data.iloc[:,:19] = data.iloc[:,:19].apply(lambda x: np.power(2,x))   # Converting log factors to linear scale
    
    #data.kcluster.unique()  # number of clusters
    #data[data['kcluster'] == 1].kcluster.unique()

    pca = None
    mu = None
    nComp = 5
    
    dataPCA = pd.DataFrame(columns = data.columns)
    dataPCAT = pd.DataFrame(columns = data.columns)
    n_clusters = 4
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
        pca = PCA(n_components = 19).fit(X)
        print(pca.components_)
        print(pca.explained_variance_)
        print(np.var(X, axis=0))
        #print(pca.singular_values_)
        print(np.dot(pca.components_[0], pca.components_[1]))
    
        # Reverse PCA
        nComp = 5
    
        samples = 1000
#        samples = 2 
        comp_1 = np.zeros(samples)
        comp_2 = np.zeros(samples) 
        m = comp_1
        for i in range(samples):
            scale_check = False
            xrand = np.zeros(shape=(1,19))
            while not scale_check:
                xrand = np.random.normal(scale=np.sqrt(pca.explained_variance_), size=(1, 19))
                Xhat = np.dot(xrand[:,:nComp], pca.components_[:nComp,:])
                Xhat += mu
                Xhat = Xhat.reshape((19,))
#                print(Xhat, x1.iloc[0,:], Xhat - x1.iloc[0,:])
                if min(Xhat) >= 0:
                    scale_check = True
                if not scale_check:
                    print("magnet factors outside bounds, re randomizing \n")        
        #    Xhat = np.zeros(19)+1
        #    Xhat = x1.iloc[0,:]
            
            resol = np.zeros(5)
            resol = cosyrun(Xhat)
#            print(new_row)
            new_row = make_db_row( Xhat, resol, cluster, data.columns )
            dataPCAtemp = pd.DataFrame(data=new_row,index=[0])
            dataPCA = dataPCA.append(dataPCAtemp,ignore_index=True)
#            new_row = make_db_row( xrand.reshape((5,)), resol, cluster, data.columns )
#            dataPCAtemp = pd.DataFrame(data=new_row,index=[0])
#            dataPCAT = dataPCAT.append(dataPCAtemp,ignore_index=True)
        plot_x, plot_y = 0, 0
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
    dataPCA.to_hdf('results_{}/best{}PCA.h5'.format(results_no,results_no),key='df')
    draw_cluster.main(data, dataPCA, 'results_{}/best{}PCA.h5'.format(results_no, results_no)) 
    dataPCA.to_hdf('../output/secar_4d_db_{}s_PCA.h5'.format(results_no,results_no),key='df')
#    print(x1.iloc[0,:],X[0,:],Xhat[0,:])
#    draw_pca(pca, x1)
    return     

def draw_pca(pca, x1):
    
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
    plt.show()
    
    # Plot PCA components # Remember to include the mean if using these components!!!
    for i in range(5):
        plt.plot(np.arange(1,20,1), pca.components_[i], label = str(i+1), linewidth=5-5*sumS[i])
    plt.legend()
    ticks = np.arange(0, len(pca.components_), 1)
    plt.xticks(ticks)
    plt.grid(True, which = 'both')
    plt.ylabel('Coefficient')
    plt.xlabel('Magnet')
    plt.show()
    
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
    plt.show()
    return


main()
