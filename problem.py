from cosy import cosyrun
from numpy import array, zeros
qNom = array([-0.39773, 0.217880+0.001472, 0.242643-0.0005+0.000729, -0.24501-0.002549, 0.1112810+0.00111, 0.181721-0.000093+0.00010-0.000096, -0.0301435+0.0001215] )

# make pygmo problem 
class optimizeRes:
    def __init__(self, dim):
        self.dim = dim

    def fitness(self, x):
        resol = -cosyrun(x)
        f = open("output1.csv","a")
        for i in range(len(x)):
            f.write("{0},".format(x[i]))
        f.write("{0}\n".format(resol))
        f.close()
        return [resol]

    def get_bounds(self):
        qLower =zeros(self.dim)
        qUpper =zeros(self.dim)
        for i in range(self.dim):
            if qNom[i] > 0:
                qLower[i] = qNom[i]*0.5
                qUpper[i] = qNom[i]*1.5
            else:
                qLower[i] = qNom[i]*1.5
                qUpper[i] = qNom[i]*0.5
        return (qLower, qUpper)

    def get_name(self):
        return "resolution optimizer"

    def get_extra_info(self):
        return "\tDimensions: " + str(self.dim)




