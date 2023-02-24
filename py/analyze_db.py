#!/usr/bin/env python3

#import commands
import pandas as pd
import numpy as np
import os,sys
import secar_utils as secar_utils

configs = secar_utils.load_configs()

# set pandas view options to print everything
#pd.set_option("max_rows", None)
#pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"

# only show best [#] since we get a lot of points
show_best = 10
kclusters = configs['clusters']
n_obj = configs['n_obj']
objectives = configs['objectives']

def main(start_i=0):

    # specify database for input
    db_out = OUTPUT_DIR + "secar_{}d_db_{}s.h5".format(n_obj, start_i)
    # initialize empty df
    df = None
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
    pareto = secar_utils.is_pareto_efficient_simple(costs)
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
    df = secar_utils.run_kmeans(df, magnet_dim, kclusters)
    df_lin = secar_utils.run_kmeans(df_lin, magnet_dim, kclusters)

    print(df.loc[df['closest']==True], df_lin.loc[df_lin['closest']==True])

    # use the replace_cluster function if want to subdivide a cluster
#    df = replace_cluster(df, 1)
    
    # only use the linear q values
    df = df_lin

    # print objective values for [show_best] number of points, sorted by FP4_BeamSpot
    show_cols = objectives
    show_cols.append('closest')   
    show_cols.append('ssobjs')
    show_cols.append('kcluster')
    print(df.loc[:show_best,show_cols].sort_values(by=['ssobjs']))

    # sort df by FP4_BeamSpot values, and reindex
    df = df.sort_values(by=objectives[-1],ignore_index=True)

    RESULTS_DIR = "results_{}/".format(start_i)
    isExist = os.path.exists(RESULTS_DIR)
    if not isExist:
        os.makedirs(RESULTS_DIR)

    # print the magnet scale factors for the best FP4_BeamSpot points
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
    (df.loc[df['closest']==True].iloc[:,:19]).round(5).to_csv(RESULTS_DIR+'magnet_factors.csv',header=write_qnames,index=False)

    # write only the magnet values and objective values to df
    df = df.drop('pareto',axis=1)
    df_out = RESULTS_DIR+"best{}.h5".format(start_i)
    df.to_hdf(df_out,key='df')

    return df_out

if __name__=='__main__':
    # input should give the batch number (i.e. the input to optimize.py)
    inputs = sys.argv
    batch = 210
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    main(batch)
