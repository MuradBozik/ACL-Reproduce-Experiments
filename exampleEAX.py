#Author: YP
#Created: 2018-08-21
#An example of how to run our scripts on the EAX data set, with comments
#explaining what each line does to serve as some documentation for how
#to use the code.

#General helper functions we use
import helper
#The main code is in ACL
import ACL as acl


#***************************************************************************
#********************** Initialize some parameters *************************
#***************************************************************************


#All work should be done in a scenario directory, with a 'phase' 
#sub-directory. In practice, you can make a new scenario directory for
#every scenario you run, if you want. However, because I found that I often
#had sevaral scenarios that share information, like their target algorithms
#and hence their parameter configuratoin space files, I further sub-divided
#my scenarios into specific "phases" of work that I am doing within that
#target algorithm's scenario.
#There are several files that my code will assume are in specific locations.
#We will point out these dependencies to you as they come up
#For example, you will need to have a file named 'params.pcs' in your
#scenario folder that uses the SMAC syntax for specifying your parameter
#configuration space for your target algorithm. 
scenario = 'scenarios/EAX'
#The phase directory must be contained within the scenario directory. 
#We will be location for files in [scenario]/[phase]/
phase = 'slices-tsp-rue-1000-3000'

#Use 'all' to process all parameters. (Note that you cannot, therefore, have a parameter of your algorithm named 'all'
parameter = 'all'
#You can also sepcify a single, case-sensitive parameter name E.g.,
#parameter = 'Npop'

#Specify the running time cutoff used in your target algorithm runs here.
cutoff = 3600

#Specify the number of (outer) bootstrap samples you want to use. We used 1001. 
numBootstrap = 1001

#Specify the number of (inner) bootstrap samples you want to use. We used 1001.
numBootstrapPerInstance = 1001

#Specify the number of independent target algorithm runs you performed per instance here
numRunsPerInstance = 10

#You can control the verbosity of the code here
#The higher the number the more verbose it will be. I generally recommend
#using verbose = 1, so that you can see how quickly progress is being made
verbose = 1

#You can choose the significance level for the bootstrap confidence intervals here. Note that we use bonferroni multiple test correction to account for the fact the each parameter response slice is made up of multiple parameter values (i.e., we will divide the number you put here by the number of parameter values used in your parameter response slices).
alpha = 0.05


#***************************************************************************
#********************* The code starts running here ************************
#***************************************************************************

#You must have a configs.csv file located in [scenario]/[phase]/configs.csv 
#that contains configuration ids, and parameter call strings for every 
#configuration. 
#See, e.g., scenarios/EAX/slices-tsp-rue-1000-3000/configs.csv

#You must also have a running time csv file containing running times and 
#miscellaneous data for each instance and seed in your instance set.
#If you don't have the seeds saved, feel free to use 0 through N, where
#N is hte number of independent runs. The miscellanous data is optional,
#you can ommit it. For legacy reasons, our code requires an integer instance
#size to be provided in this file as well. If this doesn't make sense for 
#your scenario, you can set it to any integer (e.g., -1). All that it is 
#used for in our code is to sort the instances in the outputted parameter
#response slice csv files that we used for plotting purposes. So you can
#label the size of your instances by 0,1,2,3,... etc. to control this 
#ordering, if desired. 
#The file must be formatted as followed:
#instance name, seed, instance size, running time[, miscellaneous text]
#Use 'inf' to indicate censored runs. 
#The name and location of the file must be as follows:
#[scenario]/[phase]/data/runtimes-[configuration id].csv or .log
#See scenarios/EAX/slices/tsp-rue-1000-3000/runtimes-pNPOP-0.log for an example




print("Loading the running time data...")

#Loads the data, and stores it is [scenario]/[phase]/data.pkl for later use.

#NOTE: that there must be a running time CSV file for every configuration, 
#and a configs.csv file, as specified above.

acl.loadCSVData(scenario,phase,verbose)

print("Done loading the running time data.")



print("Performing bootstrap sampling and converting to parameter response slices...")

#Next we convert the data that we just loaded from running time files into
#parameter response slices complete with bootstrap confidence intervals
#This is done for both individual instances and instance sets here.

#NOTE: We processed a large amount of data with our original scripts, and
#derived a lot more data than we ended up presenting, as much of it ended
#up not being useful. As a result, instead of passing the data around as
#inputs and outputs, we used pickle files to store and load the data 
#at the beginning and end of each function call. This allowed us to easily
#parallelize the computation of the statistics, so that we could run
#the code more quickly, especially during debugging. :)
#I'll try to clearly document here what the inputs and outputs are for each 
#function that are hidden because of this

#Hidden input: 
    #The data.pkl file created above.

#Hidden output: (Stored in pieces in [scenario]/[phase]/slice-data/)
    #responseByParam - Contains the parameter response slices for the instance sets.
    #responseByParam['parameter name'][parameter value] = [median, lower bound, upper bound]

    #perInstanceResponseByParam - Contains the parameter response slices for the individual instances, this is also a dict.
    #perInstanceResponseByParam['parameter name'][parameter value][instance # (0-n)] = [median, lower bound, upper bound]

    #bRespByPar - We needed to perform some analysis on the bootstrap samples themselves so that we could calculate confidence intervals
    #on their derived statistics (see FDC analysis). So we also save those for later use. This corresponds to the instance set response 
    #slice bootstrap sample performance statistics.
    #bRespByPar['parameter name'][bootstrap sample #][parameter value] = performance statistic (PAR10 of instance medians)

    #perInstBRespByPar['parameter name'][instance # (0-n)][bootstrap sample #][parameter value] = performance statistic (median of bootstrap sample of independent runs)

    #To make the data more easy to manipulate, we disgarded the instance name information, but it can be useful to retain this information
    #for analysis, so we also created a list of instances, and instance sizes (if applicable in your application, we were interested in analysing them 
    #as well for other threads we were pursuing, but none of this ended up making it into our main analysis), that are in the same order
    #as the instances in this analysis.

    #instances[i] = the name of the ith instance, sorted in the order used in the above datastructures
    #instanceSizes = the size of the ith instance, sorted in the roder used in the above datstructures. 

acl.convertToResponseSlices(scenario,phase,parameter,cutoff,numBootstrap,numBootstrapPerInstance,numRunsPerInstance,verbose=1,alpha=alpha)

print("Done converting to parameter response slices.")



print("Saving the parameter response slices to CSV files (e.g., for plotting)...")

#Create some CSV files for plotting/manual inspection

#Hidden input: respByParam and perInstRespByParam, from above
acl.prepPlotSliceCSVs(scenario,phase,parameter,verbose)

print("Done creating CSV files.")



#Now we're going to start doing the work of deriving statistics based
#on these parameter response slices and confidence intervals. We will
#Create a collection of statistics stored in a dict. We initialize that here
paramStats = {}

print("Calculating statistics for individual instance response slices...")

#NOTE: You don't actually have to pass in an empty dict at this point, 
#That's the default value for paramStats, But since we want to keep adding
#results to the same dict, it's best practice to initliaze one yourself,
#and then continue passing it into the following functions, they will
#return the same dict with all the information in it from before, plus
#the results of their computations.

#Hidden inputs:
    #respByParam
    #perInstRespByParam
paramStats = acl.checkPerInstanceProperties(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)

print("Done calculating statistics for individual instance response slices.")



#Next, we're going to calculate a similar set of statistics for the instance
#set parameter response slices.

print("Calculating instance set parameter response slice statistics.")

#Hidden inputs:
    #respByParam
    #perInstRespByParam

paramStats = acl.checkInstanceSetProperties(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)

print("Done calculating instance set statistics.")



#There is one additional type of analysis that we have not yet performed
#That is, FDC analysis, which requires the raw bootstrap samples, that was
#not used by the above functions. 

print("Calculating the FDC scores...")

#Hidden inputs:
    #bRespByParam
    #perInstBRespByParam

paramStats = acl.calFitnessDistances(scenario,phase,parameter,cutoff,paramStats,alpha,verbose)

print("Done calculating the FDC scores.")



#Finally, we need to save the results to a file.
#There are functions for saving the paramStats data to pkl files, and
#you can load old data while merging it with an existing paramStats
#file too. You can find these functions in ACL.py, they can be usefull
#for running things in parallel. Here, however, we will simply call
#the one functions to save paramStats to a pkl file, and to write the
#results to a csv file. If you do split things in parallel, make sure
#that you use some kind of locking mechanism so that you don't read and
#write to the files as the same time.

print("Saving the calculated statistics to a pickle file.")

acl.saveParamStats(scenario,phase,paramStats)

print("Done saving them to a pickle file.")
print("Writing the results to csv files.")

acl.writeParamStats(scenario,phase,paramStats,alpha)

print("Done writing the results to csv files.")

#If everything ran successfully, you should now find all of the computed 
#statistics and csv plot files in [scenario]/[phase]/responses/

print("The derived statistics can be found in " + scenario + '/' + phase + "responses.")

#Have fun! :)

