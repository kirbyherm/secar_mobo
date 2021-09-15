#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python
import pandas as pd
import numpy as np
import os

pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

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

PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
OUTPUT_PREFIX = OUTPUT_DIR + 'output_4f_moead_FP2_FP3_150_'
show_best = 84

db_out = OUTPUT_DIR + "secar_4d_db_90s.h5"
#db_out = OUTPUT_DIR + "better_than_nominal.h5"
df = None

objectives = ['FP2_res','FP2_e_xangle','FP3_res','FP3_e_xangle']
if os.path.exists(db_out):
    print('opening existing db')
    df = pd.read_hdf(db_out)    

print(df[objectives])

#df = df.loc[(df['FP2_res'] < max_obj) & (df['FP2_e_xangle'] < max_obj) & (df['FP3_res'] < max_obj) & (df['FP3_e_xangle'] < max_obj)]
#df = df.loc[()]

costs = df[['FP2_res','FP2_e_xangle','FP3_res','FP3_e_xangle']]
#print(costs)
costs = np.array(costs)
#print(costs)
pareto = is_pareto_efficient_simple(costs)
df['pareto'] = pareto
print(np.count_nonzero(pareto) )
df = (df.loc[(df['pareto']==True)])
df['ssobjs'] = np.sqrt(df['FP2_res']**2+df['FP2_e_xangle']**2+df['FP3_res']**2+df['FP3_e_xangle']**2)
df = df.sort_values(by='ssobjs',ignore_index=True)
df = df.loc[df['ssobjs'] < df['ssobjs'][show_best]]
df = df.sort_values(by='FP2_res',ignore_index=True)
print(df[['FP2_res','FP2_e_xangle','FP3_res','FP3_e_xangle']])
print(np.power(2,df.iloc[0,:11]))
df = df.iloc[:,:15]
df.to_hdf('test.h5',key='df')
#df.to_hdf('test.h5')#,header=False,index=False)
