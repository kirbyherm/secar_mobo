#!/usr/bin/env python3

from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

#import commands
import pandas as pd
import numpy as np
import os,sys,copy
import json
import subprocess as commands

def load_configs(configFile = 'config.json'):
    
    f = open(configFile)
    configDict = {}
    configs = json.load(f)
    keys = list(configs.keys())
    for i in range(len(keys)):
        if keys[i]=='fNominal':
            tempNom = configs['fNominal']
            plotNom = copy.deepcopy(tempNom)
            for j in tempNom.keys():
                if "_res" in j:
                    tempNom[j] = 1/tempNom[j]
#                elif j == 4:
#                    plotNom[j] = 1 * plotNom[j]
            configDict['fNominal_plot'] = plotNom
            configDict['fNominal'] = tempNom
        else:
            configDict[keys[i]] = configs[keys[i]]
    f.close()
    return configDict

configs = load_configs()

def check_fireside(): 
    #check if using fireside, i.e. whether scratch space is writeable
    computer = commands.run(['hostname'], capture_output=True).stdout.strip().decode('utf8','strict')
    if computer in ['n101','n102','n103','n104']:
        return True
    else:
        return False

# function i found which constructs a pareto front from an input set of points
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

# kmeans algorithm
def run_kmeans(df, magnet_dim=configs['magnet_dim'], clusters=configs['clusters']):

    X = df.iloc[:,:magnet_dim]
    kmeans = KMeans(n_clusters=clusters,random_state=0,n_init=10).fit(X)
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, X)
    df['kcluster'] = kmeans.labels_
    df['closest'] = False
    df.loc[closest,'closest'] = True

    return df

# expand a cluster into several more clusters
def replace_cluster(df, replace_cluster_i):

    # identify the cluster to be expanded
    replace_cluster_i = 1
    # select only these cluster points, and drop the cluster and closest id info
    df2 = df.loc[df['kcluster']==replace_cluster_i]
    df2 = df2.drop('kcluster',1) 
    df2 = df2.drop('closest',1) 
    df2 = df2.reset_index(drop=True)

    # rerun a clustering on just this cluster
    df2 = run_kmeans(df2, magnet_dim, kclusters)

    # add the results back to the original df
    df2['kcluster'] += 3
    df = df.loc[df['kcluster'] != replace_cluster_i]
    for i in range(3-replace_cluster_i):
        df['kcluster'].mask(df['kcluster'] == (replace_cluster_i+i+1),replace_cluster_i+i,inplace=True) 
    df = pd.concat([df,df2])
    df = df.reset_index(drop=True)

    return df
