#!/usr/bin/env python3

#!/mnt/simulations/secarml/soft/anaconda3/bin/python
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python


#import commands
import pandas as pd
import numpy as np
import os,sys
from sklearn.cluster import KMeans
from sklearn.metrics import pairwise_distances_argmin_min

# set pandas view options to print everything
#pd.set_option("max_rows", None)
#pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"

# only show best [#] since we get a lot of points
show_best = 10
kclusters = 4

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
def run_kmeans(df, magnet_dim=magnet_dim, clusters=clusters):

    X = df.iloc[:,:magnet_dim]
    kmeans = KMeans(n_clusters=clusters,random_state=0).fit(X)
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

def main(start_i=batch):

    # specify database for input
    db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
    # initialize empty df
    df = None
    # list the objective names
    objectives = ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    # read df
    if os.path.exists(db_out):
        print('opening existing db')
        df = pd.read_hdf(db_out)    
    
    # drop na if any
    df = df.dropna()
    
    # restrict the df to only the points that fit the problem constraints
    #   (can also change this to any value, e.g. 1 to show only better than nominal)
    max_obj = 3
    query_txt = '' 
    for i in range(len(objectives)):
        query_txt += objectives[i] + "<{}".format(max_obj)
        if i < len(objectives)-1:
            query_txt+="&"
    df = df.query(query_txt)
    
    # get costs and pass to pareto function
    costs = df[objectives]
    costs = np.array(costs)
    pareto = is_pareto_efficient_simple(costs)
    # add pareto column to df
    df['pareto'] = pareto
    print("pareto front points: {}".format(np.count_nonzero(pareto) ))
    # restrict df to only those points on the pareto front
    df = (df.loc[(df['pareto']==True)])
    df = df.reset_index(drop=True)

    # additional code which provides an alternate way of sorting/filtering df, based on sum-squares of objs
    df_objs = df[objectives]
    df_objs['ssobjs'] = (df_objs**2).sum(axis=1)
    df['ssobjs'] = df_objs['ssobjs']
    df = df.sort_values(by='ssobjs',ignore_index=True)

    quads = df.columns
    for q in range(len(quads)):
        print(q, quads[q])
        if "q" in quads[q]:
            magnet_dim = q+1

    df_lin = df.copy()
    for i in range(magnet_dim):
        df_lin.iloc[:,i] = df_lin.iloc[:,i].apply(lambda x: np.power(2,x))
    df = run_kmeans(df, magnet_dim, kclusters)
    df_lin = run_kmeans(df_lin, magnet_dim, kclusters)

    print(df.loc[df['closest']==True], df_lin.loc[df_lin['closest']==True])

    # use the replace_cluster function if want to subdivide a cluster
#    df = replace_cluster(df, 1)
    
    # only use the linear q values
    df = df_lin

    # print objective values for [show_best] number of points, sorted by FP4_BeamSpot
    print(df.loc[:show_best,['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot','closest','ssobjs','kcluster']].sort_values(by=['ssobjs']))

    # sort df by FP4_BeamSpot values, and reindex
    df = df.sort_values(by='FP4_BeamSpot',ignore_index=True)

    # print the magnet scale factors for the best FP4_BeamSpot points
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
    (df.loc[df['closest']==True].iloc[:,:19]).round(5).to_csv('magnet_factors.csv',header=write_qnames,index=False)

    # write only the magnet values and objective values to df
    df = df.drop('pareto',1)
    df.to_hdf('best{}.h5'.format(start_i),key='df')

    return

if __name__=='__main__':
    # input should give the batch number (i.e. the input to optimize.py)
    inputs = sys.argv
    batch = 210
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    main(batch)
