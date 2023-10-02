#!/usr/bin/env python3

import secar_utils as secar_utils
import make_profiles, draw_cluster_inverse
import draw_full, analyze_db, plot_tsne
import make_db

#import commands
import pandas as pd
import numpy as np
import os,sys
import subprocess as commands
import shutil
import tarfile

configs = secar_utils.load_configs()

# set pandas view options to print everything
#pd.set_option("max_rows", None)
#pd.set_option("max_columns", None)

# set important directories and files
PYGMO_DIR = "../"
OUTPUT_DIR = PYGMO_DIR + "output/"
FOX_DIR = PYGMO_DIR + 'fox/'
COSY_DIR = PYGMO_DIR + 'COSY10.0/'

# only show best [#] since we get a lot of points
kclusters = configs['clusters']
n_obj = configs['n_obj']
objectives = configs['objectives']
pca_bool = configs['pca_run']

def main(start_i=0, pca_bool=False):

    make_db.main(start_i)
    results_h5 = analyze_db.main(start_i, pca_bool)
    make_profiles.main(results_h5, pca_bool)
    draw_full.main(results_h5, start_i)
    draw_cluster_inverse.main(results_h5, results_h5) 
    plot_tsne.plot_tsne_linear(results_h5, results_h5) 
#    plot_tsne.plot_tsne(results_h5, results_h5) 

    print("\npython analysis complete, now drawing the optics profiles with cosy\n")
    os.chdir('results_{}/profiles/'.format(start_i))
    for i in range(kclusters+2):
        cmd = '../../' + COSY_DIR + 'cosy'
        cosyFilename = "pygmoCosy{}".format(i)
        foxFilename = "{}.fox".format(cosyFilename)
        output = commands.run([cmd ,foxFilename], capture_output=True)
        os.remove("{}.lis".format(cosyFilename))
        shutil.move("pic001.pdf", "X{}.pdf".format(i))
        shutil.move("pic002.pdf", "Y{}.pdf".format(i))
        cosyFilename = "pygmoCosy{}_DE".format(i)
        foxFilename = "{}.fox".format(cosyFilename)
        output = commands.run([cmd ,foxFilename], capture_output=True)
        shutil.move("pic001.pdf", "X{}_DE.pdf".format(i))
        os.remove("{}.lis".format(cosyFilename))
        os.remove("pic002.pdf")
        os.remove("RKLOG.DAT")
    os.chdir("../../")
    with tarfile.open("results_{}.tar.gz".format(start_i), "w:gz") as tar:
        filelist = [os.path.join(dp, f) for dp, dn, fn in os.walk(os.path.expanduser("results_{}".format(start_i))) for f in fn]
        for each in filelist:
            tar.add(each)

    return

if __name__=='__main__':
    # input should give the batch number (i.e. the input to optimize.py)
    inputs = sys.argv
    batch = 568
    print(inputs)
    if len(inputs) > 1:
        batch = int(inputs[1])
    if len(inputs) > 2:
        pca_bool = int(inputs[2])
    main(batch, pca_bool)

