SECAR beam tuning optimizer using pygmo

cosy.py 
    defines cosyrun(input) that runs a COSY simulation 
    given input magnet settings, 
    and returns a beam resolution value at FP1

problem.py 
    defines a pygmo UDP which calls cosyrun as its fitness

optimize.py
    defines main() which implements the pygmo archipelago
    and runs the optimization evolution
