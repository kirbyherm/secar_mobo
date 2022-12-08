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
from kmeans_cluster import run_kmeans

# set pandas view options to print everything
pd.set_option("max_rows", None)
pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"

scale_factor = np.array([0.916096997902245,1.0682652577138,0.994493371138475,0.93966084500023,1.05275223744803,1.06042964393537,1.00437784795672,0.973641800379054,1.07533403645974,1.06881462007463,1.05902890235334,1.05329541257734,0.998902975441088,1.06217562289834,1.03384085684119,1.00944081324584, 0.944682833032244,0.937425053447303,1.0784587034454])
scale_factor = np.zeros(len(scale_factor))+1

Qnom = np.array([-0.40033,0.219852,0.2552369,-0.246677876,0.11087109,0.175336731,-0.0268214976,-0.14859,0.2855,-0.0335,0.149432825,-0.182,0.1910,0.12900,-0.1380,0,0,0,0])
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

# only show best [#] since we get a lot of points
show_best = 10
batch = 210
kclusters = 4

def replace_cluster(df, replace_cluster_i):
    replace_cluster_i = 1
    df2 = df.loc[df['kcluster']==replace_cluster_i]
    df2 = df2.drop('kcluster',1) 
    df2 = df2.drop('closest',1) 
    df2 = df2.reset_index(drop=True)

    df2 = run_kmeans(df2, magnet_dim, kclusters)
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
#    if start_i == 340 or start_i == 380:
#        max_obj = 5e0
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
    if start_i in [280, 350, 360] or start_i > 379 or start_i == 0:
        for i in range(magnet_dim):
#            df.iloc[:,i] = df.iloc[:,i].apply(lambda x: np.log2(np.power(2,x)*scale_factor[i]))
            df_lin.iloc[:,i] = df_lin.iloc[:,i].apply(lambda x: np.power(2,x))
    
    df = run_kmeans(df, magnet_dim, kclusters)
    df_lin = run_kmeans(df_lin, magnet_dim, kclusters)

    print(df.loc[df['closest']==True], df_lin.loc[df_lin['closest']==True])
#    df = replace_cluster(df, 1)
    df = df_lin

    # print objective values for [show_best] number of points, sorted by FP4_BeamSpot
    print(df.loc[:show_best,['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot','closest','ssobjs','kcluster']].sort_values(by=['ssobjs']))
    # sort df by FP4_BeamSpot values, and reindex
    df = df.sort_values(by='FP4_BeamSpot',ignore_index=True)
    # print the magnet scale factors for the best FP4_BeamSpot points
    write_qnames = ['q1','q2','q3','q4','q5','q6','q7','q8','q9','q10','q11','q12','q13','q14','q15','h1','h2','h3','o1']
#    (scale_factor*np.power(2,df.loc[df['closest']==True].iloc[:,:19])).round(5).to_csv('magnet_factors.csv',header=write_qnames,index=False)
    (df.loc[df['closest']==True].iloc[:,:19]).round(5).to_csv('magnet_factors.csv',header=write_qnames,index=False)
    # write only the magnet values and objective values to df
#    print(df.columns)
    df = df.drop('pareto',1)
    #df = df.iloc[:,:19]
    df.to_hdf('best{}.h5'.format(start_i),key='df')
    return

if __name__=='__main__':
    inputs = sys.argv
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    main(batch)
