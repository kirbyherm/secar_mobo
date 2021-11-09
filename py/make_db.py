#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# import commands
import numpy as np
import pandas as pd
import os, sys

# specify important directories and names
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
generations = 750
population_size = 70
batch = 230
# specify number of magnets
magnet_dim = 6

def main(gens=generations, batch_id=210):

    OUTPUT_PREFIX = OUTPUT_DIR + 'output_4f_moead_FP2_FP3_{}_'.format(gens)
    # set up columns for dataframe
    quads = []
    for i in range(magnet_dim):
        quads.append("q{}".format(i+1))
    columns = quads
    columns.append("FP1_res")
#    columns.append("FP2_res")
#    columns.append("FP3_res")
    columns.append("MaxBeamWidth")
    columns.append("FP4_BeamSpot")
    
    # i run batches in 10s, so i specify the first id of the batch
    start_i = batch_id 
    end_i = start_i + 10
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
    better_than_nominal = 0
    for i in range(start_i, end_i):
        df_new = pd.read_csv('{}{}.csv'.format(OUTPUT_PREFIX,i),names=columns)
        if np.median(df_new.iloc[-int(len(df_new.index)*0.1):,-4:]) < 1e2:
            converged += 1
        max_obj = 1
        if len(df_new.loc[(df_new['FP1_res'] < max_obj) & (df_new['MaxBeamWidth'] < max_obj) & (df_new['FP4_BeamSpot'] < max_obj)].index) > len(df_new.index)*0.01:
#        if len(df_new.loc[(df_new['FP2_res'] < max_obj) & (df_new['FP3_res'] < max_obj) & (df_new['MaxBeamWidth'] < max_obj) & (df_new['FP4_BeamSpot'] <max_obj)].index) > len(df_new.index)*0.01:
            better_than_nominal += 1
#        print(df_new)
        # append to existing df
        if i > start_i: # or os.path.exists(db_out):
            df = df.append(df_new,ignore_index=True)
        # initialize df
        else:
            df = df_new
    
    print("percent done: {:.2f}%".format(len(df.index)/((generations+1)*population_size)*10))
    # write df to h5
    df.to_hdf(db_out,key='df')
    max_obj = 1e9
    # check for solutions strictly better than nominal (all objs < 1)
#    df = df.loc[(df['FP2_res'] < max_obj) & (df['FP3_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] <max_obj)]
    df = df.loc[(df['FP1_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df_new['FP4_BeamSpot'] < max_obj)]
    converged_points = len(df.index)
    max_obj = 1
    # check for solutions strictly better than nominal (all objs < 1)
    df = df.loc[(df['FP1_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] < max_obj)]
#    df = df.loc[(df['FP2_res'] < max_obj) & (df['FP3_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] <max_obj)]
    print("converged islands: ", converged, "\nconverged points: ", converged_points, "\nbetter than nominal islands: ", better_than_nominal, "\nbetter than nominal points: ", len(df.index))
    return

if __name__=='__main__':
    inputs = sys.argv
    print(inputs)
    if len(inputs) > 1:
        generations = int(inputs[1])
    if len(inputs) > 2:
        batch = int(inputs[2])
    main(generations,batch)
