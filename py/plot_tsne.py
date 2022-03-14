#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# Author: Narine Kokhlikyan <narine@slice.com>
# License: BSD

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from matplotlib.ticker import NullFormatter
from sklearn import manifold, datasets
from time import time

n_samples = 150
n_components = 2


def plot_tsne(filename):

    (fig, subplots) = plt.subplots(2, 4, figsize=(12, 6))
    perplexities = [5, 30, 50, 100, 300, 500, 1000]
    
    X, y = datasets.make_circles(
        n_samples=n_samples, factor=0.5, noise=0.05, random_state=0
    )
    
    df = pd.read_hdf(filename)
    print(df)
    X = np.array(df[['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','q16','q17','q18','q19']])
    y = np.array(df['kcluster'])
    z = np.array(df['closest'])
    number_of_clusters = np.max(y+1)
    colors = list(plt.get_cmap('tab20').colors)
    clusters = []
    for i in range(number_of_clusters):
        clusters.append(y == i)
    
#    red = y == 3
#    green = y == 0
#    blue = y == 2
#    yellow = y == 4
#    black = y == 1
    
    close = z == True
    
    
    print (X,y)
    
    ax = subplots[0][0]
    for i in range(number_of_clusters):
        ax.scatter(X[clusters[i],0], X[clusters[i],1], c=(colors[i]))
#    ax.scatter(X[red, 0], X[red, 1], c="c")
#    ax.scatter(X[green, 0], X[green, 1], c="m")
#    ax.scatter(X[blue, 0], X[blue, 1], c="r")
#    ax.scatter(X[yellow, 0], X[yellow, 1], c="y")
#    ax.scatter(X[black, 0], X[black, 1], c="k")
    ax.scatter(X[close,0], X[close,1], marker='x', c='b')
    ax.xaxis.set_major_formatter(NullFormatter())
    ax.yaxis.set_major_formatter(NullFormatter())
    plt.axis("tight")
    
    plot_x, plot_y = 0, 0
    for i, perplexity in enumerate(perplexities):
        plot_x += 1
        if plot_x > 3:
            plot_y+=1
            plot_x=0    
    
        ax = subplots[plot_y][plot_x]
    
        t0 = time()
        tsne = manifold.TSNE(
            n_components=n_components,
            init="random",
            random_state=0,
            perplexity=perplexity,
            learning_rate=10.0,
            n_iter=5000,
        )
        Y = tsne.fit_transform(X)
        t1 = time()
        print("circles, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))
        ax.set_title("Perplexity=%d" % perplexity)
        for j in range(number_of_clusters):
            ax.scatter(Y[clusters[j],0], Y[clusters[j],1], c=(colors[j]))
#        ax.scatter(Y[red, 0], Y[red, 1], c="c")
#        ax.scatter(Y[green, 0], Y[green, 1], c="m")
#        ax.scatter(Y[blue, 0], Y[blue, 1], c="r")
#        ax.scatter(Y[yellow, 0], Y[yellow, 1], c="y")
#        ax.scatter(Y[black, 0], Y[black, 1], c="k")
        ax.scatter(Y[close,0], Y[close,1], marker='x', c='b')
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.axis("tight")
    
    # Another example using s-curve
    #X, color = datasets.make_s_curve(n_samples, random_state=0)
    #
    #ax = subplots[1][0]
    #ax.scatter(X[:, 0], X[:, 2], c=color)
    #ax.xaxis.set_major_formatter(NullFormatter())
    #ax.yaxis.set_major_formatter(NullFormatter())
    #
    #for i, perplexity in enumerate(perplexities):
    #    ax = subplots[1][i + 1]
    #
    #    t0 = time()
    #    tsne = manifold.TSNE(
    #        n_components=n_components,
    #        init="random",
    #        random_state=0,
    #        perplexity=perplexity,
    #        n_iter=300,
    #    )
    #    Y = tsne.fit_transform(X)
    #    t1 = time()
    #    print("S-curve, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))
    #
    #    ax.set_title("Perplexity=%d" % perplexity)
    #    ax.scatter(Y[:, 0], Y[:, 1], c=color)
    #    ax.xaxis.set_major_formatter(NullFormatter())
    #    ax.yaxis.set_major_formatter(NullFormatter())
    #    ax.axis("tight")
    #
    #
    ## Another example using a 2D uniform grid
    #x = np.linspace(0, 1, int(np.sqrt(n_samples)))
    #xx, yy = np.meshgrid(x, x)
    #X = np.hstack(
    #    [
    #        xx.ravel().reshape(-1, 1),
    #        yy.ravel().reshape(-1, 1),
    #    ]
    #)
    #color = xx.ravel()
    #ax = subplots[2][0]
    #ax.scatter(X[:, 0], X[:, 1], c=color)
    #ax.xaxis.set_major_formatter(NullFormatter())
    #ax.yaxis.set_major_formatter(NullFormatter())
    #
    #for i, perplexity in enumerate(perplexities):
    #    ax = subplots[2][i + 1]
    #
    #    t0 = time()
    #    tsne = manifold.TSNE(
    #        n_components=n_components,
    #        init="random",
    #        random_state=0,
    #        perplexity=perplexity,
    #        n_iter=400,
    #    )
    #    Y = tsne.fit_transform(X)
    #    t1 = time()
    #    print("uniform grid, perplexity=%d in %.2g sec" % (perplexity, t1 - t0))
    #
    #    ax.set_title("Perplexity=%d" % perplexity)
    #    ax.scatter(Y[:, 0], Y[:, 1], c=color)
    #    ax.xaxis.set_major_formatter(NullFormatter())
    #    ax.yaxis.set_major_formatter(NullFormatter())
    #    ax.axis("tight")
    #
    #
    plt.savefig(filename+'_tsne.png')
    return


if __name__=='__main__':
    inputs = sys.argv
    plot_tsne(inputs[1])


