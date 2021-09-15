from cosy import cosyrun
from numpy import array, zeros, multiply, power

PYGMO_DIR = '../'
FOX_DIR = PYGMO_DIR + 'fox/'
OUTPUT_DIR = PYGMO_DIR + 'output/'

# make pygmo problem 
class optimizeRes:
    def __init__(self, dim, out="out"):
        self.dim = dim
        self.out = OUTPUT_DIR + out

    def fitness(self, x):
        pass_x = power(zeros(self.dim)+2.0,x)
        resol = cosyrun(pass_x)# [0:2]
        f = open(self.out,"a")
        for i in range(len(x)):
            f.write("{0},".format(x[i]))
        for i in range(len(resol)):
            if i == len(resol)-1:
                f.write("{0}\n".format(resol[i]))
            else: 
                f.write("{0},".format(resol[i]))
        f.close()
        return resol

    def get_nobj(self):
        return 4

    def get_bounds(self):
        qLower =zeros(self.dim) - 1.0
        qUpper =zeros(self.dim) + 1.0
        return (qLower, qUpper)

    def get_name(self):
        return "resolution optimizer"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)




