#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# import commands
import pandas as pd
import os

# specify important directories and names
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
generations = 800
OUTPUT_PREFIX = OUTPUT_DIR + 'output_4f_moead_FP2_FP3_{}_'.format(generations)

# specify number of magnets
magnet_dim = 15

# set up columns for dataframe
quads = []
for i in range(magnet_dim):
    quads.append("q{}".format(i+1))
columns = quads
columns.append("FP2_res")
columns.append("FP3_res")
columns.append("MaxBeamWidth")
columns.append("FP4_BeamSpot")

# i run batches in 10s, so i specify the first id of the batch
start_i = 190
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
for i in range(start_i, end_i):
    df_new = pd.read_csv('{}{}.csv'.format(OUTPUT_PREFIX,i),names=columns)
    print(df_new)
    # append to existing df
    if i > start_i: # or os.path.exists(db_out):
        df = df.append(df_new,ignore_index=True)
    # initialize df
    else:
        df = df_new

print(df)
# write df to h5
df.to_hdf(db_out,key='df')
max_obj = 1 
# check for solutions strictly better than nominal (all objs < 1)
df = df.loc[(df['FP2_res'] < max_obj) & (df['FP3_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP4_BeamSpot'] <max_obj)]
print(df)
