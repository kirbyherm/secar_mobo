#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
import pandas as pd
import os


PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
generations = 800
OUTPUT_PREFIX = OUTPUT_DIR + 'output_4f_moead_FP2_FP3_{}_'.format(generations)
magnet_dim = 15
quads = []
for i in range(magnet_dim):
    quads.append("q{}".format(i+1))
columns = quads
columns.append("FP2_res")
#columns.append("FP2_e_xangle")
columns.append("FP3_res")
#columns.append("FP3_e_xangle")
columns.append("MaxBeamWidth")
columns.append("FP4_BeamSpot")

start_i = 170
db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
df = None
#if os.path.exists(db_out):
#    print('appending')
#    df = pd.read_hdf(db_out)    
#    print(df)
end_i = start_i + 10
for i in range(start_i, end_i):
    df_new = pd.read_csv('{}{}.csv'.format(OUTPUT_PREFIX,i),names=columns)
    print(df_new)
    if i > start_i: # or os.path.exists(db_out):
        df = df.append(df_new,ignore_index=True)
    else:
        df = df_new
print(df)
df.to_hdf(db_out,key='df')
#df = df.loc[(df['FP2_res'] < 1.0) & (df['FP2_e_xangle'] < 1.0) & (df['FP3_res'] < 1.0) & (df['FP3_e_xangle'] <1.0)]
df = df.loc[(df['FP2_res'] < 1.0) & (df['FP3_res'] < 1.0) & (df['MaxBeamWidth'] < 1.0) & (df['FP4_BeamSpot'] <1.0)]
print(df)
#df.to_csv(OUTPUT_DIR + 'better_than_nominal.csv', index=False)
#df.to_hdf(OUTPUT_DIR + 'better_than_nominal.h5',key='df')
