#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# import commands
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min
import sys, os

magnet_dim = 15 
clusters = 15


def read_file( read_list, filename ):

    return

def draw_hist(df,filename):

    qNom = np.zeros(magnet_dim)
    for j in range(clusters):
        plt.cla()
        df_cluster = df.loc[df['kcluster'] == j]
        plot_x, plot_y = 0,0
        fig2, axs2 = plt.subplots(4,4)
        for i in range(magnet_dim):
            if plot_x > 3:
                plot_x, plot_y = 0, plot_y+1 
            axs2[plot_y,plot_x] = df_cluster['q{0}'.format(i+1)].plot.hist(ax=axs2[plot_y,plot_x],bins=20,range=(-1,1))
            axs2[plot_y,plot_x].axvline( x = qNom[i], ymin=0,ymax=20,color='red',linestyle='dashed')
    #        axs[plot_y,plot_x].axvline( x = max_y[i], ymin=0,ymax=20,color='green',linestyle='dashed')
            axs2[plot_y,plot_x].axes.yaxis.set_visible(False)
            axs2[plot_y,plot_x].axes.set_xlim(-1,1)
            axs2[plot_y,plot_x].set_title("q{0}".format(i+1))
            plot_x += 1
        
        fig2.delaxes(axs2[plot_y,plot_x])
        fig2.tight_layout()
        plt.savefig(filename + ".cluster_"+str(j)+".png")
        print(df.loc[(df['closest'] == True) & (df['kcluster']==j)])

    return

def run_kmeans(df, magnet_dim=magnet_dim, clusters=clusters):

#    print(magnet_dim, clusters)
    X = df.iloc[:,:magnet_dim]
    kmeans = KMeans(n_clusters=clusters,random_state=0).fit(X)
    closest, _ = pairwise_distances_argmin_min(kmeans.cluster_centers_, X)
    df['kcluster'] = kmeans.labels_
    df['closest'] = False
    df.loc[closest,'closest'] = True
    return df

def main():

    script, filename = sys.argv
    df = pd.read_hdf(filename)
    df = run_kmeans(df)
    draw_hist(df,filename)
    print("hello world")

if __name__ == "__main__":
    main()
