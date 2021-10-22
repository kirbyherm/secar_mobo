#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python

# make sure above path points to the version of python where you have pygmo installed 
# nscl servers
#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
# hpcc servers
#!/mnt/home/herman67/anaconda3/envs/pygmo/bin/python


#import commands
import pandas as pd
import numpy as np
import os
from kmeans_cluster import run_kmeans

# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"

Qnom = np.array([-0.40033,0.219852,0.2552369,-0.246677876,0.11087109,0.175336731,-0.0268214976,-0.14859,0.2855,-0.0335,0.149432825,-0.182,0.1910,0.12900,-0.1380])
# Q1H:=0.003703;
# H1:=0.0103564;
# H2:=0.0052735{*0.5};
# H3:=-0.008774463{*1.5};
# O1:=0.031283{*2.0};

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

# only show best 100 since we get a lot of points
show_best = 100 

start_i = 180
# specify database for input
db_out = OUTPUT_DIR + "secar_4d_db_{}s.h5".format(start_i)
# initialize empty df
df = None
# list the objective names
objectives = ['FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
# read df
if os.path.exists(db_out):
    print('opening existing db')
    df = pd.read_hdf(db_out)    

# drop na if any
df = df.dropna()

# restrict the df to only the points that fit the problem constraints
#   (can also change this to any value, e.g. 1 to show only better than nominal)
max_obj = 1e9
df = df.loc[(df['FP2_res'] < max_obj) & (df['MaxBeamWidth'] < max_obj) & (df['FP3_res'] < max_obj) & (df['FP4_BeamSpot'] < max_obj)]

# get costs and pass to pareto function
costs = df[['FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']]
costs = np.array(costs)
pareto = is_pareto_efficient_simple(costs)
# add pareto column to df
df['pareto'] = pareto
print(np.count_nonzero(pareto) )
# restrict df to only those points on the pareto front
df = (df.loc[(df['pareto']==True)])
df = df.reset_index(drop=True)
print(df.iloc[:100,:4])
# additional code which provides an alternate way of sorting/filtering df, based on sum-squares of objs
#df['ssobjs'] = np.sqrt(df['FP2_res']**2+df['FP3_res']**2+df['MaxBeamWidth']**2+df['FP4_BeamSpot']**2)
#df = df.sort_values(by='ssobjs',ignore_index=True)
#df = df.loc[df['ssobjs'] < df['ssobjs'][show_best]]
df = run_kmeans(df, 15, 10)

# sort df by FP4_BeamSpot values, and reindex
df = df.sort_values(by='FP4_BeamSpot',ignore_index=True)
# print objective values for [show_best] number of points, sorted by FP4_BeamSpot
#print(df.iloc[:,15:])
# print the magnet scale factors for the best FP4_BeamSpot points
(np.power(2,df.loc[df['closest']==True].iloc[:,:15])).round(5).to_csv('magnet_factors.csv',index=False)
(Qnom * np.power(2,df.loc[df['closest']==True].iloc[:,:15])).round(5).to_csv('magnet_values.csv',index=False)
# write only the magnet values and objective values to df
print(df.columns)
df = df.drop('pareto',1)
#df = df.iloc[:,:19]
df.to_hdf('best.h5',key='df')

