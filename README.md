SECAR quad tuning optimizer using pygmo

Requirements:
    pygmo
    pandas
    numpy
    matplotlib
    (a variety of other py packages are imported, but most should be preinstalled)

Files:

    cosy.py 
        defines cosyrun(input) that runs a COSY simulation given input magnet settings 
            and returns a set of objective values
        assumes cosy is installed and can be run with:
            cosy foxFile.fox
    
    problem.py 
        defines a pygmo User-Defined Problem which calls cosyrun in its fitness evaluation 
        *pygmo_testing* uses a range of -0.1 to 0.1 in order to generate results within the bounds
    
    optimize.py
        defines main() which implements the pygmo archipelago and runs the optimization evolution
            (note that the archipelago can be a single island) 
    
    make_db.py
        defines a method for constructing a pandas hdf5 database from the csv results
        db to be made is specified within the code 
    
    view_db.py
        defines a method to read the hdf5 db and write only values which meet a criteria to a separate db
        db to be viewed is specified within the code    

    draw.py
        defines a method for plotting the results
        requires an argument of the file to be analyzed 
        e.g. run with:
            ./draw.py ../output/outputFileName.csv    

Notes:
    with proper updates to the shebang lines all scripts can simply be called with ./[script name]
    
