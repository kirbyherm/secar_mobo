#!/usr/bin/env python3

# import commands
from cosy import cosyrun
from numpy import array, zeros, multiply, power
import secar_utils as secar_utils 

configs = secar_utils.load_configs()
fox_name = configs['fox_name']
fNom = configs['fNominal']
nobj = configs['n_obj']

# set important directories
PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
OUTPUT_DIR = PYGMO_DIR + 'output/'
SCRATCH_DIR = configs['scratch_dir']

# make pygmo problem 
class optimizeRes:
    # define the output file to write solutions
    #   this is how we see all the attempted solutions
    #   and check how the algorithm performs
    # output file is defined when initiating the problem in optimize.py
    def __init__(self, dim, out="out"):
        self.dim = dim
        self.out = out

    # fitness function that is called by each iteration
    #   x is between -1 and 1, but cosyrun expects a factor
    #   between 2^-1 and 2^1
    def fitness(self, x):
        # convert to scale factor
        pass_x = power(zeros(self.dim)+2.0,x)
        # run cosy
        resol = cosyrun(pass_x, fNom, fox_name)
        # return objective values to evolve
        return resol

    # define number of objectives
    def get_nobj(self):
        return nobj

    # define bounds of x values
    #   i.e. from -1 to 1 (powers of 2)
    def get_bounds(self):
        qLower =zeros(self.dim) - 2.
        qUpper =zeros(self.dim) + 2.
        qLower[9] *= 1.5 
        qUpper[9] *= 1.5 
        return (qLower, qUpper)

    # define problem name
    def get_name(self):
        return "resolution optimizer"

    # print dimensions 
    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)




