#!/mnt/misc/sw/x86_64/all/anaconda/python3.7/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from scipy.spatial.distance import directed_hausdorff

from itertools import combinations

def calc_hausdorff(filename,filename_compare):
    
    df = pd.read_hdf(filename)
    df_compare = pd.read_hdf(filename_compare)
#    print(df)
    objectives = ['FP1_res','FP2_res','FP3_res','MaxBeamWidth','FP4_BeamSpot']
    objective_compare = 'FP4_BeamSpot'
    magnets = []
    for i in range(1,20):
        magnets.append("q{}".format(i))

#    u = np.array(df[magnets])
#    v = np.array(df_compare[magnets])
    fig, axs = plt.subplots(2,2)
    plot_y, plot_x = 0, 0
    for i in range(0,4):
        
#        print(df[objective_compare], df[objectives[i]])
        axs[plot_y, plot_x].plot(df[objective_compare], df[objectives[i]],'o',linestyle='None',markersize=0.3)
        axs[plot_y, plot_x].plot(df_compare[objective_compare], df_compare[objectives[i]],'o',linestyle='None',markersize=0.3)
        axs[plot_y, plot_x].set_ylabel(objectives[i])
        if plot_y == 1:
            axs[plot_y, plot_x].set_xlabel(objective_compare)
        u = np.array(df[[objectives[i],objective_compare]])
        v = np.array(df_compare[[objectives[i],objective_compare]])
#        axs[i].set_yscale('log')
#        axs[i].set_xscale('log')
        dh_uv = directed_hausdorff(u,v)
        dh_vu = directed_hausdorff(v,u)
        print(u.shape, v.shape, dh_uv, dh_vu)
        axs[plot_y, plot_x].plot([u[dh_uv[1]][1],v[dh_uv[2]][1]],[u[dh_uv[1]][0],v[dh_uv[2]][0]], 'x', linestyle='None',fillstyle='none')
        axs[plot_y, plot_x].plot([u[dh_vu[2]][1],v[dh_vu[1]][1]],[u[dh_vu[2]][0],v[dh_vu[1]][0]], 's', linestyle='None',fillstyle='none')
        print(objectives[i], dh_uv[0], dh_vu[0])
        plot_x += 1
        if plot_x > 1:
            plot_y += 1
            plot_x = 0
    plt.show()
    return

def toy_model():

    xdim = 10000
    u = np.zeros((int(2*xdim),2))
    v = np.zeros((int(2*xdim),2))

    r1, r2 = 1, 3
    for i in range(xdim):
        x = np.random.uniform(-r1,r1)
        u[i,0] = x
        u[i+xdim,0] = x
        u[i,1] = np.sqrt(r1**2- x**2)
        u[i+xdim,1] = -np.sqrt(r1**2-x**2)
        x *= r2/r1
        v[i,0] = x
        v[i+xdim,0] = x
        v[i,1] = np.sqrt(r2**2- x**2)
        v[i+xdim,1] = -np.sqrt(r2**2-x**2)
    print('expected difference in radii = {}'.format(abs(r1-r2)))
    gen_hausdorff = (max(directed_hausdorff(u,v)[0],directed_hausdorff(v,u)[0]))
    print('calculated hausdorff distance = {0:.2f}'.format(gen_hausdorff))

    plt.plot(u[:,0],u[:,1],marker='o',linestyle='None')
    plt.plot(v[:,0],v[:,1],marker='o',linestyle='None')
    plt.show()

    return


if __name__=='__main__':
#    toy_model()


#    inputs = sys.argv
    test_list = ['280','330','340','350']
    name_list = []
    for i in test_list:
        name_list.append('results_{0}/best{0}.h5'.format(i))
    print(name_list)
    compares = list(combinations(name_list,2))
    compares = [compares[1]]
    print(compares)

    for i in compares:
#$        print(i[0], i[1])
        calc_hausdorff(i[0],i[1])

