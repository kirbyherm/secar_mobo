#!/usr/bin/env python3

# import commands
import numpy as np
import pandas as pd
import os, sys

# specify important directories and names
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
generations = 750
population_size = 1001
batch = 230
# specify number of magnets
magnet_dim = 19

def cut_data(df, objectives, max_obj):
    if type(max_obj) in [int,float]: max_obj = np.zeros(len(objectives))+max_obj
    query_txt = '' 
    for i in range(len(objectives)):
        query_txt += objectives[i] + "<{}".format(max_obj[i])
        if i < len(objectives)-1:
            query_txt+="&"
    df_return = df.query(query_txt)
    return df_return

def main(batch_id=210, gens=generations ):

    OUTPUT_PREFIX = OUTPUT_DIR + 'output'
    if batch_id < 100 and batch_id > 0:
        OUTPUT_PREFIX = OUTPUT_DIR + 'output_4f_nsga2_FP2_FP3_{}_'.format(gens)
    # set up columns for dataframe
    objectives = ['FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    quads = []
    for i in range(magnet_dim):
        quads.append("q{}".format(i+1))
    columns = quads
    for i in range(len(objectives)):
        columns.append(objectives[i])
    
    # i run batches in 10s, so i specify the first id of the batch
    start_i = batch_id 
    end_i = start_i + 1
    if start_i > 360 and start_i < 390:
        end_i = start_i + 1
    elif start_i == 0:
        start_i = 280
        end_i = 440
    # name the output
    db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
    
    # check if we want to append to an existing df
    df = None
    #if os.path.exists(db_out):
    #    print('appending')
    #    df = pd.read_hdf(db_out)    
    #    print(df)
    
    # load results from all islands
    converged = 0
    converged_islands = []
    better_than_nominal = 0
    better_than_nominal_islands = []
    for i in range(start_i, end_i):
        if start_i == 280 and end_i == 440:
            if i > 289 and i < 410:
                continue
        df_new = pd.read_csv('{}{}.csv'.format(OUTPUT_PREFIX,i),names=columns)
        if np.median(df_new.iloc[-int(len(df_new.index)*0.1):,-len(objectives):]) < 1e2:
            converged += 1
            converged_islands.append(i)
        max_obj = 1.01
        if len(cut_data(df_new,objectives,max_obj).index) > 0:
            better_than_nominal += len(cut_data(df_new,objectives,max_obj).index)
            better_than_nominal_islands.append(i)
        # append to existing df
        if i > start_i: # or os.path.exists(db_out):
            df = df.append(df_new,ignore_index=True)
        # initialize df
        else:
            df = df_new
    
    print("completed sims: {:}\npercent done: {:.2f}%".format(len(df.index),len(df.index)/((generations+1)*population_size*(end_i - start_i))*100))
    # write df to h5
#    if start_i in [280, 350, 360]:
#        for i in range(magnet_dim):
##            df.iloc[:,i] = df.iloc[:,i].apply(lambda x: np.log2(np.power(2,x)*scale_factor[i]))
#            df.iloc[:,i] = df.iloc[:,i].apply(lambda x: np.power(2,x))
#    df.loc[:,'MaxBeamWidth'] = df.loc[:,'MaxBeamWidth'].apply(lambda x: np.power(x, 1/4.0))
    df.to_hdf(db_out,key='df')
    max_obj = 1e9
    # check for solutions strictly better than nominal (all objs < 1)
    df = cut_data(df, objectives, max_obj)
#    df = df.loc[(df['FP1_res'] < max_obj) & (df['FP2_res'] < max_obj) & (df['FP3_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] <max_obj)]
    converged_points = len(df.index)
    max_obj = 1.01
    # check for solutions strictly better than nominal (all objs < 1)
#    df = df.loc[(df['FP1_res'] < max_obj) & (df['FP2_res'] < max_obj) & (df['FP3_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] <max_obj)]
    df = cut_data(df, objectives, max_obj)
    print("converged islands: ", converged_islands, "\nconverged points: ", converged_points, "\nbetter than nominal islands: ", better_than_nominal_islands, "\nbetter than nominal points: ", len(df.index))
    return

if __name__=='__main__':
    inputs = sys.argv
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    main(batch)
