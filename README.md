# SECAR quad tuning optimizer using pygmo

## Requirements:

pygmo  
pandas  
numpy  
matplotlib  
scikit-learn  
(a variety of other py packages are imported, but most should be preinstalled)
    
the python packages can be easily installed with conda using: 
`conda env create -f environment.yml`  
then activate with:  
`conda activate secar_moead`

Additionally, change the INCLUDE path in the header of the fox files in the `./fox/` directory to point to the location of COSY.bin

## Files:

### COSY10.0/ :
    
directory holding all of the necessary scripts and binaries for executing cosy 

if cosy is already installed on your system you can ignore this and simply use your version
presuming you update the fox files to use the correct `INCLUDE` path and change `COSY_DIR` in python scripts

### fox/ :

contains all the fox files necessary for running the optimization

files should resemble `SECAR_pg_Optics.fox`, `SECAR_pg_Optics_DE.fox`,  
or `SECAR_an_Optics.fox, SECAR_an_Optics_DE.fox`
`SECAR_pg_Optics.fox` specifies the optics of just the recoils through SECAR
and prints out the dimensions at each element in addition to the dimensions for the objectives
in this case that is the X width at each FP1, 2, 3, and the beam spot size on the DSSD
`SECAR_pg_Optics_DE.fox` specifies the optics of the unreacted beam (or whatever contaminant you'd like to separate) through SECAR
    and prints out the location of the contaminant for each objective
        in this case the separation at FP1, 2, 3 

change `fox_name` in `config.json` to match the base name of your desired fox file (i.e. ``SECAR_pg_Optics``)

as above, be sure to update the INCLUDE path in the header of the fox files to point to the location of `COSY.bin`

most importantly, any new fox files need to print out X and Y dimensions at every element of the beamline
see `cosy.py` header, specifically the variables `magnet_names` and `magnet_dims` for the elements, or else use `SECAR_pg_Optics.fox` as a guide

both the files must also print out the desired objective values at their locations, for `SECAR_pg_Optics.fox` and `SECAR_pg_Optics_DE.fox`
that means the recoil widths and contaminant separations at FP1, 2, 3 and the recoil beam spot size at the DSSD

the `MaxBeamWidth` is then calculated by taking the largest width of the recoils relative to the width of each element
see `cosy.py` for specifics

### py/ : 

#### `analyze_db.py`

defines a method to read the hdf5 db and write only values which are on the pareto front and are less than `max_obj` to a separate db
also runs a kmeans-clustering algorithm on the db to generate clusters of points  

#### `config.json`

specifies most of the variables needed for running the full process
comments don't exist in json so additional keys like `fNominalan` and `fNominalpg` are unused 
they are simply there to hold the values 

#### `cosy.py`

defines cosyrun(input) that runs a COSY simulation given input magnet settings 
and returns a set of objective values
assumes cosy is installed and can be run with:
`cosy foxFile.fox`
all the interfacing between cosy and python is done here, so be sure to check that, given a specified reaction type and `fNominal` values, 
`./cosy.py` returns values which are == 1 or only off by a precision < 1e-7 or so


#### `optimize.py`
defines main() which implements the pygmo archipelago and runs the optimization evolution
(note that the archipelago can be a single island) 

#### `pca.py`
runs a pca on the specified results and then generates a number of random points using those components

#### `pipeline.py`
loads in all the necessary python scripts and runs them to generate and plot results

#### `problem.py`
defines a pygmo User-Defined Problem which calls cosyrun in its fitness evaluation 

#### `utils.py`
defines various useful functions, most notably providing the wrapper to read in the `config.json` file    

### Drawing scripts :

Several different `draw*.py` scripts exist, in various states of clarity.
I will continue to try to clean and comment these but they are less important than the analysis scripts

#### `draw.py`
defines a method for plotting the results
requires an argument of the file to be analyzed 
e.g. run with:  
`./draw.py ../output/outputFileName.csv`


### sh/ :
    
all the bash scripts for submitting jobs to fireside are here

#### `analysis_job.sh`
launches a job for `pca.py`, generating the pca points

#### slurmfiles/ :
try to write out the sbatch logs to this file so they can be easily stored and cleaned as necessary

#### `submit_job.sh`
launches a job for `optimize.py`
        

## Deprecated scripts:    

results are no longer stored to csv, and so making them into a pandas df is trivial

#### `make_db.py`  
defines a method for constructing a pandas hdf5 database from the csv results
db to be made is specified within the code 

#### `pipeline.sh`  
basic pipelining done in bash, updated to the `pipeline.py` script

## Notes:

assuming python3 is installed with conda and  
`which python3`  
points to the correct version all scripts can simply be made executable and called with:  
`./[script name] [args] `
    
