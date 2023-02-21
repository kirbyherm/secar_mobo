#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from matplotlib.ticker import NullFormatter
from sklearn import manifold, datasets
from time import time

import matplotlib

# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

colormap = matplotlib.cm.get_cmap('PuOr')

n_samples = 150
n_components = 2


def plot_corr(filename):
    
    df = pd.read_hdf(filename)
#    df = df.loc[df['kcluster']==3]
    print(df)
    magnets = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','q16','q17','q18','q19']
    magnets = magnets[:15]
    magnet_scale_factors = []
    df = df[magnets].apply(lambda x: np.power(2,x))
    for i in range(len(magnets)):
        if magnets[i] == 'q10':
            magnet_scale_factors.append([-3,3])
        else: 
            magnet_scale_factors.append([-2,2])
    plot_combos = [[0,1],[2,3],[3,4],[2,4],[4,6],[6,9],[9,14],[12,14]]

    X = np.array(df[magnets])

    
    print(df[magnets])
    print(df[magnets].corr() )
    fig, ax = plt.subplots()
    correlations = df[magnets].corr()
    for i in range(correlations.shape[0]):
        for j in range(correlations.shape[1]):
            if abs(correlations.iloc[i,j]) < 0.5:
                correlations.iloc[i,j]=0.0
        
    plt.pcolormesh(correlations,cmap=colormap)
    labels = magnets 
    ax.set_xticks(np.arange(len(magnets))+0.5)
    ax.set_xticklabels(np.arange(len(magnets))+1)
    ax.set_yticks(np.arange(len(magnets))+0.5)
    ax.set_yticklabels(np.arange(len(magnets))+1)
    ax.set_ylabel('q')
    ax.set_xlabel('q')
    ax.set_title('correlation matrix for q1-q15')
    plt.colorbar()
    plt.savefig(filename+'_corr.png')

    return


if __name__=='__main__':
    inputs = sys.argv
    plot_corr(inputs[1])


