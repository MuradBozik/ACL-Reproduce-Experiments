#Author: YP
#Created: 2018-08-21
#An example of how to run our scripts on all of the data sets.
#See exampleEAX.py for more detailed comments on what the code
#is doing.

#NOTE: Running this code will take several hours.

#General helper functions we use
import helper
#The main code is in ACL
import ACL as acl

#***************************************************************************
#********************** Initialize some parameters *************************
#***************************************************************************

#Use 'all' to process all parameters. (Note that you cannot, therefore, have a parameter of your algorithm named 'all'
parameter = 'all'
#You can also sepcify a single, case-sensitive parameter name E.g.,
#parameter = 'Npop'

#Specify the number of (outer) bootstrap samples you want to use. We used 1001. 
numBootstrap = 1001

#Specify the number of (inner) bootstrap samples you want to use. We used 1001.
numBootstrapPerInstance = 1001

#Specify the number of independent target algorithm runs you performed per instance here
numRunsPerInstance = 10

#You can control the verbosity of the code here
#The higher the number the more verbose it will be. I generally recommend
#using verbose = 1, so that you can see how quickly progress is being made
verbose = 0

#You can choose the significance level for the bootstrap confidence intervals here. Note that we use bonferroni multiple test correction to account for the fact the each parameter response slice is made up of multiple parameter values (i.e., we will divide the number you put here by the number of parameter values used in your parameter response slices).
alpha = 0.05


#***************************************************************************
#********************* The code starts running here ************************
#***************************************************************************


tsp = {'slices-tsp-rue-1000-3000':3600}
sat = {'slices-circuit-fuzz':600,'slices-BMC08':600}
mip = {'slices-cls':10000,'slices-regions200':10000}

scenarios = {'scenarios/LKH':tsp,'scenarios/EAX':tsp,'scenarios/CaDiCaL':sat,'scenarios/lingeling':sat,'scenarios/cryptominisat':sat,'scenarios/CPLEX':mip}


for scenario in scenarios.keys():
  for phase in scenarios[scenario].keys():
    cutoff = scenarios[scenario][phase]

    print('*'*50)
    print("Starting " + scenario + '/' + phase)
    print('*'*50)

    print("Loading the running time data...")
    acl.loadCSVData(scenario,phase,verbose)
    print("Done loading the running time data.")

    print("Performing bootstrap sampling and converting to parameter response slices...")
    acl.convertToResponseSlices(scenario,phase,parameter,cutoff,numBootstrap,numBootstrapPerInstance,numRunsPerInstance,verbose=1,alpha=alpha)
    print("Done converting to parameter response slices.")

    print("Saving the parameter response slices to CSV files (e.g., for plotting)...")
    acl.prepPlotSliceCSVs(scenario,phase,parameter,verbose)
    print("Done creating CSV files.")

    paramStats = {}
    print("Calculating statistics for individual instance response slices...")
    paramStats = acl.checkPerInstanceProperties(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)
    print("Done calculating statistics for individual instance response slices.")

    print("Calculating instance set parameter response slice statistics.")
    paramStats = acl.checkInstanceSetProperties(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)
    print("Done calculating instance set statistics.")

    print("Calculating the FDC scores...")
    paramStats = acl.calFitnessDistances(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)
    print("Done calculating the FDC scores.")

    print("Saving the calculated statistics to a pickle file.")
    acl.saveParamStats(scenario,phase,paramStats)
    print("Done saving them to a pickle file.")

    print("Writing the results to a csv file.")
    acl.writeParamStats(scenario,phase,paramStats,alpha)
    print("Done writing the results to a csv file.")

#If everything ran successfully, you should now find all of the computed 
#statistics and csv plot files in [scenario]/[phase]/responses/
#for all of the scenarios we studied.

#Have fun! :)

