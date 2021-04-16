from cosy import cosyrun
from numpy import array, zeros, multiply, power
qNom = array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )

output_file = "output_4f_moead_n2.csv"
# make pygmo problem 
class optimizeRes:
    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        pass_x = power(zeros(self.dim)+2.0,x)
        resol = cosyrun(multiply(pass_x,qNom))# [0:2]
        f = open(output_file,"a")
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
#        for i in range(self.dim):
#            if qNom[i] > 0:
#                qLower[i] = 0.5 # qNom[i]*0.5
#                qUpper[i] = 1.5 # qNom[i]*1.5
#            else:
#                qLower[i] = 1.5 # qNom[i]*1.5
#                qUpper[i] = 0.5 # qNom[i]*0.5
        return (qLower, qUpper)

    def get_name(self):
        return "resolution optimizer"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)




