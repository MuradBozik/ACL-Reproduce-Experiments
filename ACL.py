#Author: YP
#Created: 2017-12
#A collection of functions for reproducing the results of the 2018 PPSN paper: 
#Algorithm Configuration Landscapes: More Benign than Expected?


import numpy as np
import glob
import copy
import random 
import math
import datetime
import scipy
from scipy.spatial import ConvexHull

import helper



def Exp(x,p=1.5):
    return p**x - 1 


def realRange(centre,minVal,maxVal,n,balance=0.75,f=Exp,p=1.5,verbose=0):
    #Author: YP
    #Created: 2017-12
    #We want to find a set of logarithmically spaced numbers that, 
    #when reflected about 0 and added to the centre, contain exactly n 
    #values between minVal and maxVal and which have at most balance% of
    #the numbers on either size of 0.

    if(verbose > 1):
        print("centre: " + str(centre))
        print("min: " + str(minVal))
        print("max: " + str(maxVal))

    if(maxVal == centre or minVal == centre):
        if(verbose > 1):
            print("Setting balance = 1 (that is, removing the constraint) since the centre is at the boundary of the range")
        balance = 1


    shifts = [f(i,p) for i in range(1, int(n*balance)+1)]
    shifts2 = [-f(i,p) for i in range(1,int(n*balance)+1)]
    shifts.extend(shifts2)
    shifts.append(0)
    shifts = np.array(sorted(shifts))
    
    maxDist = max(abs(maxVal-centre),abs(minVal-centre))
    minDist = min(abs(maxVal-centre),abs(minVal-centre))
    
    maxAlpha = maxDist/(f(int(n/2.0*balance)+1,p))
    minAlpha = minDist/(f(int(n*balance)+1,p))

    contained = -1
    while(not (contained == n)):
        #Pick a random weight, alpha, for the shifts 
        alpha = minAlpha + (maxAlpha-minAlpha)*random.random()
        values = centre + alpha*shifts
        inRange = np.logical_and(minVal < values, values < maxVal)
        contained = np.sum(inRange)
        #Update the min and max ranges.
        if(contained < n): 
            maxAlpha = alpha
        elif(contained > n): 
            minAlpha = alpha

        if(verbose > 0): 
            print('-'*60)
            print('Tried alpha: ' + str(alpha))
            print('Tried values: ' + str(values[np.where(inRange)]))
            print('Number of values above: ' + str(np.sum(values[np.where(inRange)] > centre)))
            print('Number of values below: ' + str(np.sum(values[np.where(inRange)] < centre)))
            print('Number of values: ' + str(contained))

            with open('values.log','w') as f_out:
                for val in values[np.where(inRange)]:
                    f_out.write(str(val) + '\n')

    lastVal = -float('inf')
    for val in values[np.where(inRange)]:
        if(lastVal == val):
            print('[Warning]: Underflow has occured!')
        lastVal = val

    values = values[np.where(inRange)]
    #The centre value may have changed a bit due to floating point errors. Reset it to the exact value
    ind = np.argmin((values - centre)**2)
    values = list(values)
    values[ind] = centre
       
    return values


def integerRange(centre,minVal,maxVal,n,balance=0.75,f=Exp,p=1.5,verbose=0):
    #Author: YP
    #Created: 2017-12

    if(maxVal-minVal+1 <= n):
        return range(minVal,maxVal+1)

    #print(3)

    if(verbose > 1):
        print(max(min(centre,maxVal-0.49),minVal+0.49))

    values = realRange( max(min(centre,maxVal-0.49),minVal+0.49) ,minVal,maxVal,n,balance,f,p,verbose)

    #print(values)

    newValues = []
    random.shuffle(values)
    for val in sorted(values,key=lambda k: (k-centre)**2):
        if(verbose > 1):
            print(val)
        newValues, newVal = addValue(val,newValues,minVal,maxVal,centre)
        if(verbose > 1):
            print(str(val) + ' -> ' + str(newVal))

    values = sorted(values)
    newValues = sorted(newValues)

    if(verbose > 0):
        with open('values.log','w') as f_out:
            for i in range(0,len(values)):
                if(len(newValues) > i):
                    f_out.write(str(values[i]) + ', ' + str(newValues[i]) + '\n')
                else:
                    f_out.write(str(values[i]) + ', NaN\n')


    return newValues



def addValue(val,values,minVal,maxVal,centre):
    #Author: YP
    #Created: 2017-12
    #A helper function for the realRange and integerRange functions
    if(len(values) == maxVal-minVal+1):
        #Every value has already been added.
        return values
    elif(int(round(val)) not in values):
        #add the value
        values.append(int(round(val)))
        return values, int(round(val))
    #print('-- val=' + str(val) + '; centre=' + str(centre))
    if(val > centre):
        #We already have this value, so we increment the value by 1
        #if we can and try that one. However, if we have already reached the
        #boundary then we try adding the nearest value below the centre.
        if(int(round(val)) < maxVal):
            return addValue(int(round(val+1)),values,minVal,maxVal,centre)
        else:
            return addValue(centre-1,values,minVal,maxVal,centre)
    else:
        #We are below the centre and we already have this value, so we
        #Try decrementing the value by 1 and adding that value. If we hit
        #the boundary, we go the other way.
        if(int(round(val)) > minVal):
            return addValue(int(round(val-1)),values,minVal,maxVal,centre)
        else:
            return addValue(centre+1,values,minVal,maxVal,centre)


def convertToResponseSlices(scenario,phase,parameter,cutoff,numBootstrap,numBootstrapPerInstance,numRunsPerInstance,verbose=0,alpha=0.05):
    #Author: YP
    #last updated: 2018-08-21
    #I used this to collect reformat the data into a more useful form, i.e., from the raw running time data I collected
    #I used this to calculate the bootstrap sample confidence intervals and medians, both for the individual instance
    #response slices and the instance set response slices. 
    #This function is a bit of a beast, sorry..
    #This function stores the data in the new format for you so that you don't need to keep re-calling it every time.
    #This can save you some time, since this function isn't exactly fast.

    #The stored output:
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


    #helper.mkdir(scenario + '/' + phase + '/responses')
    if(verbose > 0):
        print(scenario + '/' + phase) 

    data = helper.loadObj(scenario + '/' + phase,"data")
    #data['configuration ID'][instance size]['instance name'][seed]

    if(verbose > 0):
        print("Loaded the data")

    splitData = {}

    #Build a mapping from configuration ID to parameter values
    configs = getConfigs(scenario,phase)
    values = {}


    for configId in configs.keys():
        param = '-'.join(configId.split('-')[:-1])[1:].lower()
        configString = ' ' + configs[configId]
        for p in configString.split(' -'):
            p = p.strip()
            if(p.split(' ')[0].lower() == param):
                values[configId] = float(p.split(' ')[1].strip("'"))

    counts = {}

    for configId in configs.keys():
       for size in data[configId].keys():
           if(size not in counts.keys()):
               counts[size] = {}
           for inst in data[configId][size].keys():
               if(inst not in counts[size].keys()):
                   counts[size][inst] = 0
               if(len(data[configId][size][inst]) == numRunsPerInstance):
                   counts[size][inst] += 1
 
    for size in counts.keys():
        for inst in counts[size].keys():
            if(counts[size][inst] == len(configs.keys())):
                continue
            for configId in configs.keys():
                if(size in data[configId].keys()):
                    if(inst in data[configId][size].keys()):
                        #Remove this instance because data for it has not been collected for every configuration.
                        del data[configId][size][inst]

    if(verbose > 0):
        print("Removed instances that do not have runs for every configuration")

    for configId in configs.keys():
        param = '-'.join(configId[1:].split('-')[:-1]).lower()
        if(param not in splitData.keys()):
            splitData[param] = {}
   
        #print(values.keys())
        #print(data.keys())
        splitData[param][values[configId]] = data[configId]

    if(verbose > 0):
        print("split the data by parameter")

    if(parameter == 'all'):
        parameter = list(splitData.keys())
    else:
        parameter = [parameter.lower()]

    #splitData[param][value][instance size]['instance name'][seed]
    responseByParam = {}
    perInstanceResponseByParam = {}
    bRespByPar = {} #bRespByPar[param][k][val] = the PAR10 of the kth bootstrap sample for the value val of parameter param.
    perInstBRespByPar = {} #perInstBRespByPar[param][i][k][val] = the median running time of the kth bootstrap sample on the ith instance for the value val of the parameter param
    bMeansByValue = {}

    #We save the bootstrap indices for each bootstrap sample so that each configuration is applied to the same bootstrap sample.
    savedIndsMajor = []
    savedIndsMinor = []
    savedIndices = []
    

    for param in parameter:
        #Flatten all the data for this parameter
        flatData = {}
        values = sorted(splitData[param].keys())
        for value in values:
            flatData[value] = []
            for size in sorted(splitData[param][value].keys()):
                for inst in sorted(splitData[param][value][size].keys()):
                    instData = []
                    for seed in splitData[param][value][size][inst].keys():
                        instData.append(min(splitData[param][value][size][inst][seed][0],cutoff))
                    flatData[value].append(instData)  
            flatData[value] = np.array(flatData[value])    
   
        if(verbose > 0):
            print("Flattened the data for parameter " + param)

        #----------------------------------------------------
        #We have loaded and flatten the raw data into a more convenient format to work with 
        #flatData[value][i][s] = running time of the sth run on the ith instance when the parameter is set to value         
        #Now we can begin performing the bootstrap analysis to create the confidence intervals    

        #Sorry about the confusing names...
        response = {}
        perInstanceResponses = {}
        bmeansByValue = []
        perInstanceBstats = {} 
        perInstBRespByPar[param] = {}
        bRespByPar[param] = {}

        #Apply Bonferroni multiple test correction
        alphaBC = alpha/len(values)

        print("Values: " + str(values))

        for value in values:
            if(verbose > 0):
                print("Working on value " + str(value) + "...")
                #print(len(flatData[value]))
            perInstanceResponses[value] = {}

            databPI = [] #databPI[i][k] = the median running time for the kth bootstrap sample of the ith instance running times on this value.
            for i in range(0,len(flatData[value])):

                #if(verbose > 1):
                #    print("Instance: " + str(i))
                #    print("instData: " + str(flatData[value][i]))
                #Get the running times for this instance
                instData = flatData[value][i]
                if(i >= len(savedIndices)):
                    #Create a set of bootstrap sample indices
                    savedIndices.append(np.random.randint(0,len(instData),(numBootstrapPerInstance,len(instData))))
                indices = savedIndices[i]
                #sample the data using the indices
                bInstData = np.array(instData)[indices]
                #Get the per-instance medians for each sample
                bInstMedians = np.apply_along_axis(lambda sample: helper.calStatistic(sample,'median'),1,bInstData)               
                #Add the per instance bootstrap sample medians to the bootstrap sample data array.
                databPI.append(bInstMedians)

                #Calculate a confidence interval
                perInstanceResponses[value][i] = [helper.calStatistic(bInstMedians, 'median'), helper.calStatistic(bInstMedians, 'q' + str(alphaBC/2)), helper.calStatistic(bInstMedians, 'q' + str(1-alphaBC/2))]

                if(perInstanceResponses[value][i][0] > perInstanceResponses[value][i][2]):
                    print("It really is wrong here...")

                if(i not in perInstanceBstats.keys()):
                    perInstanceBstats[i] = []
                perInstanceBstats[i].append(bInstMedians)

                if(i not in perInstBRespByPar[param].keys()):
                    perInstBRespByPar[param][i] = {}
                for k in range(0,len(bInstMedians)):
                    if(k not in perInstBRespByPar[param][i].keys()):
                        perInstBRespByPar[param][i][k] = {}
                    perInstBRespByPar[param][i][k][value] = bInstMedians[k]

       
      
            #print("Created per instance bootstrap samples")
        
            databPI = np.array(databPI)
            #databPI[i][j] = the jth bootstrap sample median running time for the ith instance
            bdata = []
            #For each outer bootstrap sample
            for sample in range(0,numBootstrap):
                #We want to always use the same bootstrap sample of the indices for every 
                #configuration, so that the differences we observe between parameter values 
                #for a single bootstrap sample set are not simply due to the variance of 
                #having taken two different bootstrap samples.
                if(sample >= len(savedIndsMajor)):
                    #sample instance indices with replacement
                    savedIndsMajor.append(np.random.randint(0,len(databPI),(len(databPI))))
                    #sample per instance bootstrap sample median indices with replacement
                    savedIndsMinor.append(np.random.randint(0,numBootstrapPerInstance,(len(databPI))))

                indsMajor = savedIndsMajor[sample]
                #sample per instance bootstrap sample median indices with replacement
                indsMinor = savedIndsMinor[sample]
                #Create the bootstrap sample and add it to the list of samples
                bdata.append(databPI[(indsMajor,indsMinor)])

                #if(sample%10 == 0):
                #    print("Created " + str(sample) + " bootstrap samples")

            #calculate the list of bootstrap sample means
            bmeans = np.apply_along_axis(lambda sample: helper.calStatistic(sample, 'PAR10', cutoff),1,bdata)
            bmeansByValue.append(bmeans)

            med = helper.calStatistic(bmeans, 'median')
            lo =  helper.calStatistic(bmeans, 'q' + str(alphaBC/2))
            up = helper.calStatistic(bmeans, 'q' + str(1-alphaBC/2))
 
            response[value] = [med,lo,up]

            for k in range(0,len(bmeans)):
                if(k not in bRespByPar[param].keys()):
                    bRespByPar[param][k] = {}
                bRespByPar[param][k][value] = bmeans[k]

            #databPIbyValue.append(databPI)

        responseByParam[param] = response
        perInstanceResponseByParam[param] = perInstanceResponses

    if(verbose > 0):
        print("Sorting instances for tracking purposes...")

    instances = []
    instanceSizes = []

    values = sorted(splitData[parameter[0]].keys())
    for value in values:
        for size in sorted(splitData[parameter[0]][value].keys()):
            instances.extend(sorted(splitData[parameter[0]][value][size].keys()))
            instanceSizes.extend([size]*len(splitData[parameter[0]][value][size].keys()))

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

    if(verbose > 0):
        print("Saving the results...")

    #Save the data for future recall   
    saveSliceData(scenario, phase, responseByParam, perInstanceResponseByParam, bRespByPar, perInstBRespByPar, instances, instanceSizes)



def saveSliceData(scenario, phase, responseByParam, perInstanceResponseByParam, bRespByPar, perInstBRespByPar, instances, instanceSizes):
    #Author: YP
    #Last updated: 2018-08-20
    #I used this to store the reformated data, processed into bootstrap confidence intervals and probabilities

    storeDir = scenario + '/' + phase + '/slice-data'
    helper.mkdir(storeDir)

    with helper.acquireLock(storeDir + '/slice-lock.log'):
        for param in responseByParam.keys():
            helper.saveObj(storeDir, responseByParam[param], 'aggregate-response-' + param)
            helper.saveObj(storeDir, perInstanceResponseByParam[param], 'per-instance-responses-' + param)
            helper.saveObj(storeDir, bRespByPar[param], 'bootstrap-aggregate-response-' + param)
            helper.saveObj(storeDir, perInstBRespByPar[param], 'bootstrap-per-instance-responses-' + param)

        helper.saveObj(storeDir, instances, 'instances')
        helper.saveObj(storeDir, instanceSizes, 'instance-sizes')


def loadSliceData(scenario, phase, param):
    #Author: YP
    #last updated: 2018-08-20
    #A short-cut created to load the reformatted data rather than
    #reformatting it every time.  

    perInstanceResponseByParam = {}
    responseByParam = {}

    storeDir = scenario +'/' + phase + '/slice-data/'

    if(param == 'all'):
        iFiles = glob.glob(storeDir + '/per-instance-responses*.pkl')
        aFiles = glob.glob(storeDir + '/aggregate-response*.pkl')
    else:
        iFiles = [storeDir + '/per-instance-responses-' + param + '.pkl']
        aFiles = [storeDir + '/aggregate-response-' + param + '.pkl']

    for f in iFiles:
        param = '-'.join(f.split('/')[-1].split('-')[3:])[:-4]
        perInstanceResponseByParam[param] = helper.loadObj(storeDir, 'per-instance-responses-' + param)

    for f in aFiles:
        param = '-'.join(f.split('/')[-1].split('-')[2:])[:-4]
        responseByParam[param] = helper.loadObj(storeDir, 'aggregate-response-' + param)

    #instances = helper.loadObj(storeDir, 'instances')
    #instanceSizes = helper.loadObj(storeDir, 'instance-sizes')

    return responseByParam, perInstanceResponseByParam


def loadBootstrapSliceData(scenario, phase, param):
    #I'm copying the above function, so the variable names are missleading here.
    
    perInstanceResponseByParam = {}
    responseByParam = {}

    storeDir = scenario +'/' + phase + '/slice-data/'

    if(param == 'all'):
        iFiles = glob.glob(storeDir + '/bootstrap-per-instance-responses*.pkl')
        aFiles = glob.glob(storeDir + '/bootstrap-aggregate-response*.pkl')
    else:
        iFiles = [storeDir + '/bootstrap-per-instance-responses-' + param + '.pkl']
        aFiles = [storeDir + '/bootstrap-aggregate-response-' + param + '.pkl']

    for f in iFiles:
        param = '-'.join(f.split('/')[-1].split('-')[4:])[:-4]
        perInstanceResponseByParam[param] = helper.loadObj(storeDir, 'bootstrap-per-instance-responses-' + param)

    for f in aFiles:
        param = '-'.join(f.split('/')[-1].split('-')[3:])[:-4]
        responseByParam[param] = helper.loadObj(storeDir, 'bootstrap-aggregate-response-' + param)

    #instances = helper.loadObj(storeDir, 'instances')
    #instanceSizes = helper.loadObj(storeDir, 'instance-sizes')


    return responseByParam, perInstanceResponseByParam


def loadInstanceData(scenario,phase):
    #Author: YP
    #Created: 2018-08-21

    storeDir = scenario +'/' + phase + '/slice-data/'

    instances = helper.loadObj(storeDir, 'instances')
    instanceSizes = helper.loadObj(storeDir, 'instance-sizes')

    return instances, instanceSizes


def prepPlotSliceCSVs(scenario,phase,parameter,verbose=0):
    #Author: YP
    #Last updated: 2018-08-21
    #I believe this was used to create the CSV files I originally
    #was using to plot the parameter response slices.  However, I think
    #I also ended up not using them for any of my final plots. But,
    #I'll need to check on that. 


    helper.mkdir(scenario + '/' + phase + '/responses')

    responsesByParam, perInstanceResponseByParam =  loadSliceData(scenario, phase, parameter)

    for param in responsesByParam.keys():
        values = sorted(responsesByParam[param].keys())
        perInstanceResponses = perInstanceResponseByParam[param]
        response = responsesByParam[param]

        with open(scenario + '/' + phase + '/responses/per-instance-responses-' + param + '.csv','w') as f_out:
            f_out.write("#Parameter value, instance 1 median of median running time, instance 1 95% lower bound on median running time, instance 1 95% upper bound on median running time, instance 2 median of median running time, ... \n")
            for value in values:
                f_out.write(str(value))
                for i in range(0,len(perInstanceResponses[value])):
                    f_out.write(', ' + ', '.join([str(perInstanceResponses[value][i][j]) for j in range(0,3)]))
                f_out.write('\n')

        with open(scenario + '/' + phase + '/responses/response-' + param + '.csv','w') as f_out:
            f_out.write("#Parameter value, median of Mean running time, 95% lower bound on mean running time, 95% upper bound on mean running time\n")
            for value in values:
                f_out.write(str(value) + ', ' + ', '.join([str(response[value][i]) for i in range(0,3)]) + '\n')




def checkInstanceSetProperties(scenario,phase,parameter,cutoff,paramStats={},alpha=0.05,verbose=0):
    #Author: YP
    #Last udpated: 2018-08-20
    #This function calculates some of the derived statistics of the instance set parameter response slices.

    respByParam, PerInstRespByParam = loadSliceData(scenario, phase, parameter)

    unimodality = {} 
    convexity = {} 
    trivialConvexity = {}
    nearFlatness = {}
   
    for param in respByParam.keys():
        if(param not in paramStats.keys()):
            paramStats[param] = {}

        values = sorted(respByParam[param].keys())
        if(verbose > 0):
            print("Working on: " + str(param))

        response = []
        for v in values:
            #NOTE: I am reordering the data here for convenience in the functions below.
            response.append([respByParam[param][v][1], respByParam[param][v][0], respByParam[param][v][2]])

        #Check to see if the response is unimodal
        bitonic = checkModality(response,verbose) == 1
        #Check to see if the response is convex
        convex, triviallyConvex = checkConvexity(values,response,cutoff*10,verbose)
        #Check to see if the response is nearly flat, i.e., not interesting
        nearlyFlat = checkFlat(values,response,cutoff,verbose)
       
        unimodality[param] = bitonic
        convexity[param] = convex
        trivialConvexity[param] = triviallyConvex
        nearFlatness[param] = nearlyFlat

    for param in sorted(respByParam.keys()):
        paramStats[param]['Unimodal'] = unimodality[param]
        paramStats[param]['Convex'] = convexity[param]
        paramStats[param]['Trivially Convex'] = trivialConvexity[param]
        paramStats[param]['Interesting'] = not nearFlatness[param]
              
    return paramStats



def checkPerInstanceProperties(scenario,phase,parameter,cutoff,paramStats={},alpha=0.05,verbose=0):
    #Author: YP
    #Last updated: 2018-08-20
    #This function derives data from the reformatted, bootstrap interval-based parameter response slices. 

    responsesByParam, perInstanceResponsesByParam = loadSliceData(scenario, phase, parameter)

    unimodality = {}
    convexity = {}
    trivialConvexity = {}
    numModes = {}
    nearFlatness = {}

    for param in perInstanceResponsesByParam.keys():
        if(param not in paramStats.keys()):
            paramStats[param] = {}
        values = sorted(perInstanceResponsesByParam[param].keys())
        pir = perInstanceResponsesByParam[param]

        unimodal = 0
        convexCount = 0
        triviallyConvexCount = 0
        totalCount = 0
        numModes[param] = []
        
        unimodality[param] = []
        convexity[param] = []
        trivialConvexity[param] = []
        nearFlatness[param] = []

        if(verbose > 1):
            print("*"*60 + '\n' + str(param))
        
        for inst in range(0,len(pir[values[0]])):
            response = []
            for v in values:
                #NOTE: I am reordering the data from (mid, low, up) to (low, mid, up) here because it is
                #a more convenient storage format for this function....
                response.append([pir[v][inst][1], pir[v][inst][0], pir[v][inst][2]])
 
            if(verbose > 1):
                print("instance: " + str(inst))

            #Check to see if it is convex (and whether or not it is trivially so, because there are too many censored datapoints to perform a real test) 
            convex, triviallyConvex = checkConvexity(values,response,cutoff,verbose)
            #Check to see if it is flat, i.e., un-interesting
            nearlyFlat = checkFlat(values,response,cutoff,verbose) 
            #Count the number of modes, and get their locations for later use
            modeCount = checkModality(response,verbose)
          
            if(modeCount == 1):
                 unimodal += 1
            if(convex):
                 convexCount += 1
            if(triviallyConvex):
                 triviallyConvexCount += 1
            totalCount += 1
            unimodality[param].append(modeCount==1)
            numModes[param].append(modeCount)
            convexity[param].append(convex)
            trivialConvexity[param].append(triviallyConvex)
            nearFlatness[param].append(nearlyFlat)

        if(verbose > 0):
            print("Parameter " + param + " had " + str(unimodal) + '/' + str(totalCount) + ' unimodal instance responses and ' + str(convexCount) + '/' + str(totalCount) + ' convex instance responses.') 

        if("Instance Stats" not in paramStats[param].keys()):
            paramStats[param]['Instance Stats'] = {}
        for j in range(0,len(unimodality[param])):
            if(j not in paramStats[param]['Instance Stats'].keys()):
                paramStats[param]['Instance Stats'][j] = {}
            paramStats[param]['Instance Stats'][j]['Unimodal'] = unimodality[param][j]
            paramStats[param]['Instance Stats'][j]['Convex'] = convexity[param][j]
            paramStats[param]['Instance Stats'][j]['Trivially Convex'] = trivialConvexity[param][j]
            paramStats[param]['Instance Stats'][j]['Interesting'] = not nearFlatness[param][j]
            paramStats[param]['Instance Stats'][j]['Number of Modes'] = numModes[param][j]
           
    for param in sorted(unimodality.keys()):
        paramStats[param]['% Instances Unimodal'] = float(sum(unimodality[param]))/len(unimodality[param])*100
        paramStats[param]['% Instances Convex'] = float(sum(convexity[param]))/len(convexity[param])*100
        paramStats[param]['% Instances Trivially Convex'] = float(sum(trivialConvexity[param]))/len(convexity[param])*100
        paramStats[param]['% Instances Interesting'] = 100 - float(sum(nearFlatness[param]))/len(nearFlatness[param])*100
        paramStats[param]['Average number of modes per instance'] = float(sum(numModes[param]))/len(numModes[param])
                  
    return paramStats



def checkModality(response,verbose=0):
    #Author: Yasha Pushak
    #First Created: March 1st, 2018
    #Checks the k-tonic piecewise linear line with minimal k that fits
    #within the confidence interval.
    #response = [[low0, mid0, up0], [low1, mid1, up1], ..., [lowN, midN, upN]]
    ##NOTE: The ordering that I am using (low, mid, up) is different from what I used before (mid, low, up).
    phase = -1 #downward
    phaseSwitches = 0
    pt = response[0][2]

    path = [] 
    pathPhases = []

    for i in range(0,len(response)):
        if( (pt - response[i][1-phase])*phase <= 0):
            #The next bound is still moving in the right direction
            #update the point to the next bound's value
            pt = response[i][1-phase]
        elif( (pt - response[i][1+phase])*phase <= 0):
            #The next bound is less constrained in this direction than the
            #last; however, we can draw a flat line and still be contained
            #within this interval, so we do not need to change phases
            pass
        else: 
            #We are on the downward/upward pass, but 
            #the next lower/upper bound is higher/lower than our current position,
            #so we have to switch to the upward/downward phase now
            phase *= -1
            pt = response[i][1-phase]
            phaseSwitches += 1
        path.append(pt)
        pathPhases.append(phase)

    numModes = int(phaseSwitches/2) + 1

    if(verbose > 1):
        print("It has " + str(numModes) + " modes.")
        if(verbose > 3):
            print("The path taken was: " + str(path))

    #Find the argminimums of the modes
    best = float('inf')
    argbest = -1
    curPhase = -1
    modeLocations = []
    for i in range(0,len(pathPhases)):
        if(pathPhases[i] == -1):
            if(response[i][1] <= best):
                best = response[i][1]
                argbest = i
            curPhase = -1
        elif(curPhase == -1):
            curPhase = 1
            modeLocations.append(argbest)
            argbest = -1
            best = float('inf')
    if(curPhase == -1):
        modeLocations.append(argbest)

    return numModes


def checkFlat(values,response,cutoff,verbose=0):
    #Author: Yasha Pushak
    #Created: 2018-03-21
    #Originally this checked to see if you can draw a straight line through the confidence intervals
    #However, we found that this straight line test eliminated some response slices
    #that were visually interesting. As a result, we created a "nearly-flat" test
    #to see if the overlap between their intervals is small enough. See the paper
    #for details on this test (Note that in the paper it described in terms of log
    #running times, here we use the expononent to equivallently perform the test 
    #using the regular coordinate space.)

    maxMinRes = -float('inf')
    minMaxRes = float('inf')

    intSizes = []
    for res in response:
        maxMinRes = max(maxMinRes,res[0])
        minMaxRes = min(minMaxRes,res[2])
        intSizes.append(res[2]/res[0])

    meanIntSize = sum(intSizes)/len(intSizes)
    minOverlapIntSize = minMaxRes/maxMinRes    
    
    if(verbose > 1):
        print("Mean interval size: " + str(meanIntSize))
        print("Minimum interval overlap size: " + str(minOverlapIntSize))
        if(minOverlapIntSize <= meanIntSize**0.5):
            print("It is not nearly flat, i.e., it is interesting")
        else:
            print("It is nearly flat, i.e., it is not interesting")

    return minOverlapIntSize >= meanIntSize**0.5
  

def checkConvexity(values,response,cutoff,verbose=0):
    #Author: Yasha Pushak
    #First Created: February 7th, 2018
    #Checks to see if the convex hull of the upper bounds ever falls below the lower bounds
    #if so, then this response is not convex.
    #The last boolean returned indicates whether or not the respones is trivially convex becasue everything is censored. 
    # 
    #response = [[low0, mid0, up0], [low1, mid1, up1], ..., [lowN, midN, upN]]
    ##NOTE: The ordering that I am using (low, mid, up) is different from what I used before (mid, low, up).
    #values must be in sorted order.

    if(not values == sorted(values)):
        raise Exception("Values not in sorted order.")

    if(verbose > 3):
        print(response) 
    #print(cutoff)

    upperBound = []
    for i in range(0,len(values)):
        if(response[i][2] < cutoff):
            upperBound.append([values[i],response[i][2]])

    if(len(upperBound) <= 2):
        if(verbose > 3):
            print("It is trivially convex because we only have " + str(len(upperBound)) + " non-censored upper bounds.")
        return True, True

    #Force the smallest and largest parameter values to be a part of the convex hull
    upperBound.append([upperBound[0][0],cutoff])
    upperBound.append([upperBound[-2][0],cutoff])

    if(verbose > 3):
        print("Upperbound: " + str(upperBound))

    upperBound = np.array(upperBound)
    hull = list(upperBound[ConvexHull(upperBound).vertices])
    
    if(verbose > 3):
        print("Hull: " + str(list(hull)))
 
    strippedHull = [[-float('inf'),float('inf')]]
    #remove the extra points
    for i in range(0,len(hull)):
        if(hull[i][1] < cutoff):
            strippedHull.append(hull[i])
    strippedHull.append([float('inf'),float('inf')])
    hull = strippedHull

    #print(hull)
    #Sort by the value
    hull = sorted(hull,key=lambda k: k[0])

    convex = True

    #For each point i in the lower bound check to see if it is above the convex hull
    j = 1
    for i in range(0,len(values)):
        if(response[i][0] >= cutoff):
            #Skip this parameter value since it is censored.
            continue
        #Find the first point in the hull that has a value greater than or equal to point i
        while(values[i] > hull[j][0]):
            j += 1
        
        if(values[i] == hull[j][0]):
            #We just need to check the height of the hull here since they are at the same location
            if(response[i][0] > hull[j][1]):
                convex = False
                if(verbose > 1):
                    print("It is not convex.")
                    print("Certificate: " + str([values[i], response[i][0]]) + " > " + str(list(hull[j])))
                break
        else:
            #We need to interpolate between this and the previous point on the convex hull.
            #We use linear interpolation. 
            hx1 = hull[j-1][0]; hy1 = hull[j-1][1]
            hx2 = hull[j][0]; hy2 = hull[j][1]
            lbx = values[i]; lby = response[i][0]

            hy = hy1 + (hy2-hy1)*((lbx-hx1)/(hx2-hx1))

            if(response[i][0] > hy):
                convex = False
                if(verbose > 1):
                    print("It is not convex.")
                    print("Certificate: " + str([values[i], response[i][0]]) + " > " + str([values[i],hy]))
                break

    if(convex and verbose > 1):
        print("It is convex.")

    #print(hull)

    return convex, False


def calFitnessDistances(scenario, phase, parameter, cutoff, paramStats={},alpha=0.05,verbose=0):
    #Author: Yasha Pushak
    #Created: 2018-03-09?
    #Calculates the median and 95% confidence interval (over bootstrap samples)
    #for the fitness-distance correlation coefficients for the parameter responses, 
    #and the average instance parameter responses.   

    if(verbose > 1):
        print("Working on calFitnessDistnaces")
 
    bRespByParam, perInstBRespByParam = loadBootstrapSliceData(scenario, phase, parameter)
    #bRespByParam[param][k][val] = the PAR10 of the kth bootstrap sample for the value val of parameter param.
    #perInstBRespByParam[param][i][k][val] = the median running time of the kth bootstrap sample on the ith instance for the value val of the parameter param

    for param in bRespByParam.keys():
        if(verbose > 1):
            print("Working on " + param)
        if(param not in paramStats.keys()):
            paramStats[param] = {}

        rFDC = []
        for k in bRespByParam[param].keys(): #Iterate over the bootstrap samples
            r = calFitDistCor(bRespByParam[param][k],float('inf'))
            if(r is not None):
                rFDC.append(r)

        paramStats[param]['Fitness Distance Correlation (Median)'] = helper.calStatistic(rFDC,'median')
        paramStats[param]['Fitness Distance Correlation (Q' + str(100*alpha/2) + ')'] = helper.calStatistic(rFDC,'q' + str(alpha/2))
        paramStats[param]['Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')'] = helper.calStatistic(rFDC,'q' + str(1-alpha/2))

        rFDCi = []
        
        for i in sorted(perInstBRespByParam[param].keys()):
            if(verbose > 1):
                print("Working on instance " + str(i))
            rFDCi.append([])
            for k in perInstBRespByParam[param][i].keys():
                r = calFitDistCor(perInstBRespByParam[param][i][k],cutoff)
                if(r is not None):
                    rFDCi[i].append(r)

        meanrFDCi = []
        for k in perInstBRespByParam[param][0].keys():
            if(verbose > 0):
                print("Working on bootstrap sample " + str(k))
            tmp = []
            for i in sorted(perInstBRespByParam[param].keys()):
                if(k < len(rFDCi[i])): #We skipped over bootstrap samples where the instance response was censored for all, or all but one value.
                    tmp.append(rFDCi[i][k])
            meanrFDCi.append(helper.calStatistic(tmp,'mean'))

        #print(len(rFDCi))
        #print(len(rFDCi[0]))
        #rFDCi = np.array(rFDCi)

        #meanrFDCi = np.apply_along_axis(lambda k: helper.calStatistic(k,'mean'),0,rFDCi)
             
        #print("Length of meanrFDCi")
        #print(len(list(meanrFDCi)))

        paramStats[param]['Average Instance Fitness Distance Correlation (Median)'] = helper.calStatistic(meanrFDCi,'median')
        paramStats[param]['Average Instance Fitness Distance Correlation (Q' + str(100*alpha/2) + ')'] = helper.calStatistic(meanrFDCi,'q' + str(alpha/2))
        paramStats[param]['Average Instance Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')'] = helper.calStatistic(meanrFDCi,'q' + str(1-alpha/2)) 

        if('Instance Stats' not in paramStats[param].keys()):
            paramStats[param]['Instance Stats'] = {} 

        if(verbose > 1):
            print("Calculating final statistics")
        for i in sorted(perInstBRespByParam[param].keys()):
            if(i not in paramStats[param]['Instance Stats'].keys()):
                paramStats[param]['Instance Stats'][i] = {}
            if(len(rFDCi[i]) == 0):
                #print(rFDCi[i])
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Median)'] =  float('NaN')
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Q' + str(100*alpha/2) + ')'] =  float('NaN')
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')'] =  float('NaN')
            else:
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Median)'] =  helper.calStatistic(rFDCi[i],'median')
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Q' + str(100*alpha/2) + ')'] =  helper.calStatistic(rFDCi[i],'q' + str(alpha/2))
                paramStats[param]['Instance Stats'][i]['Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')'] =  helper.calStatistic(rFDCi[i],'q' + str(1-alpha/2))



        #print(paramStats[param].keys())


    return paramStats
        

def calFitDistCor(response,cutoff):
    #Author: Yasha Pushak
    #Created 2018-03-13
    #Calculates the fitness distance correlation as described in SLS Foundations and Applications
    #where the fitness of a solution is it's (PAR10) running time divided by the best known solution's running time
    #response[value] = time

    best = float('inf')
    argBest = -1
    values = list(response.keys())
    random.shuffle(values) #break ties uniformly at random
    countComplete = 0
    for val in values: 
        if(response[val] >= cutoff):       
            response[val] = cutoff*10
        if(not response[val] == cutoff*10):
            countComplete += 1
        if(response[val] <= best):
            best = response[val]
            argBest = val

    if(countComplete <= 1):
        #We don't have enough data to perform this analysis. Skip this instance.
        return None
    else:
        allSame = True
        sameVal = response[values[0]]
        for val in values:
            if(not sameVal == response[val]):
                allSame = False
        if(allSame):
            #All of the running times are exactly the same. This happens because of machine precise errors for measuring running times. 
            #We are going to treat this the same way as a censored instance and skip it.
            return None



    d = []
    g = []

    for val in values:
        d.append(abs(val - argBest))
        g.append(response[val]/best)

    m = len(d)
    d = np.array(d)
    g = np.array(g)

    cov = np.dot((g - np.mean(g)),(d - np.mean(d)))/(m-1)
    sigg = np.sqrt(np.dot(g - np.mean(g), g - np.mean(g))/(m - 1))
    sigd = np.sqrt(np.dot(d - np.mean(d), d - np.mean(d))/(m - 1))

    if(math.isnan(cov/(sigg*sigd))):
        print("*"*60)
        print("fitness distance correlation is nan for the following:")
        print("response: " + str(response))
        print("Covariance: " + str(cov))
        print("sigg: " + str(sigg))
        print("sigd: " + str(sigd))
        print("g: " + str(g))
        print("d: " + str(d)) 

    return cov/(sigg*sigd)
    


def mergeAndWriteResults(scenario,phase,paramStats,logFile,alpha):
    
    if(logFile is not None):
        f_out = open(logFile,'a',0)
        f_out.write('-'*10 + datetime.datetime.now().strftime("%I:%M%p on %B %d, %Y") + '-'*10 + '\n')
        f_out.write(scenario + '/' + phase + '/Completed collecting results.\n')
    else:
        f_out = None


    with helper.acquireLock(scenario + '/' + phase + '/lock.log'):
        paramStats = loadAndMerge(scenario,phase,paramStats,f_out)
        saveParamStats(scenario,phase,paramStats,f_out)
        writeParamStats(scenario,phase,paramStats,alpha)

    if(f_out is not None):
        f_out.write("Finished writing to the Parameter-Statistics*.csv files\n")
        f_out.close()




def saveParamStats(scenario,phase,paramStats,f_out=None):
    helper.saveObj(scenario + '/' + phase + '/slice-data/',paramStats, 'Parameter-Statistics')
    if(f_out is not None):
        f_out.write("Saved updated results in pickle file\n")



def loadAndMerge(scenario,phase,newParamStats,f_out=None):

    if(helper.isFile(scenario + '/' + phase + '/slice-data/Parameter-Statistics.pkl')):
        paramStats = helper.loadObj(scenario + '/' + phase + '/slice-data/','Parameter-Statistics')
        if paramStats is None:
            return newParamStats
    else:
        return newParamStats

    #print(paramStats)

    for param in newParamStats.keys():
        if(param not in paramStats.keys()):
            paramStats[param] = {}
        for stat in newParamStats[param].keys():
            if(not stat == 'Instance Stats'):
                if(stat in paramStats.keys() and f_out is not None):
                    f_out.write("Writing over results for " + stat + " for " + param + '\n')
                #print(paramStats.keys())
                #print(newParamStats.keys())
                paramStats[param][stat] = newParamStats[param][stat]
            else:
                if('Instance Stats' not in paramStats[param].keys()):
                    paramStats[param]['Instance Stats'] = {}
                for inst in newParamStats[param]['Instance Stats'].keys():
                    if(inst not in paramStats[param]['Instance Stats'].keys()):
                        paramStats[param]['Instance Stats'][inst] = {}
                    for stat in newParamStats[param]['Instance Stats'][inst].keys():
                        if(stat in paramStats[param]['Instance Stats'][inst].keys() and f_out is not None and int(inst) == 1):
                            f_out.write("Writing over instance statistics for " + stat + " for " + param + '\n')
                        paramStats[param]['Instance Stats'][inst][stat] = newParamStats[param]['Instance Stats'][inst][stat]

    if(f_out is not None):
        f_out.write("Finished loading and merging results in pickle file\n")

    return paramStats
                    


def writeParamStats(scenario,phase,paramStats,alpha):

    instances, instanceSizes = loadInstanceData(scenario,phase)

    helper.mkdir(scenario + '/' + phase + '/responses') 
 
    stats = ['Interesting','Unimodal','Convex','Trivially Convex','% Instances Interesting','% Instances Unimodal','% Instances Convex','% Instances Trivially Convex','Average number of modes per instance','Fitness Distance Correlation (Median)','Fitness Distance Correlation (Q' + str(100*alpha/2) + ')','Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')','Average Instance Fitness Distance Correlation (Median)','Average Instance Fitness Distance Correlation (Q' + str(100*alpha/2) + ')','Average Instance Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')']

    with open(scenario + '/' + phase + '/responses/Parameter-Statistics.csv','w') as f_out:
        f_out.write("#[1] Parameter")
        i = 2
        for stat in stats:
            f_out.write(', [' + str(i) + '] ' + str(stat))
            i += 1
        f_out.write('\n')
        for param in sorted(paramStats.keys()):
            #print(paramStats[param].keys())
            f_out.write(param)
            for stat in stats:
                if(stat not in paramStats[param].keys()):
                    f_out.write(', Incomplete...')
                elif(type(paramStats[param][stat]) in [bool, np.bool_]):
                    if(paramStats[param][stat]):
                        f_out.write(', 1')
                    else:
                        f_out.write(', 0')
                else:
                    f_out.write(', ' + str(paramStats[param][stat]))
            f_out.write('\n')

    stats = ['Interesting','Unimodal','Convex','Trivially Convex','Number of Modes','Fitness Distance Correlation (Median)','Fitness Distance Correlation (Q' + str(100*alpha/2) + ')','Fitness Distance Correlation (Q' + str(100*(1-alpha/2)) + ')']



    for param in paramStats.keys():
        if('Instance Stats' not in paramStats[param].keys()):
            continue
        with open(scenario + '/' + phase + '/responses/Parameter-Statistics-' + param + '.csv','w') as f_out:
            f_out.write('#[1] Instance #, [2] Instance name')
            i = 3
            for stat in stats:
                f_out.write(', [' + str(i) + '] ' + str(stat))
                i += 1
            f_out.write('\n')

            for inst in range(0,len(paramStats[param]['Instance Stats'])):
                f_out.write(str(inst) + ', ' + instances[inst])
                for stat in stats:
                    if(inst not in paramStats[param]['Instance Stats'].keys() or stat not in paramStats[param]['Instance Stats'][inst].keys()):
                        f_out.write(', Incomplete...')
                    else:
                        f_out.write(', ' + str(paramStats[param]['Instance Stats'][inst][stat]))
                f_out.write('\n')



def getConfigs(scenario,phase):
    #Author: Yasha Pushak
    #Last updated: February 20th, 2017
    #A helper function that parses the configs.csv files to create a dict that
    #maps configuration ID to parameter call string.

    configs = {} 
    with open(scenario + '/' + phase + '/configs.csv') as f_in:
        for line in f_in:
            if('#' in line[0]):
                continue
            pId = line.split(',')[0].strip()
            params = line.split(',')[1].strip()
            configs[pId] = params

    return configs



def loadCSVData(scenario,phase,verbose=1):
    #Author: YP
    #Created: 2018-21-08
    #A helper function that loads CSV file running time data
    #into a format that is used as input for the rest of this code.

    configs = getConfigs(scenario,phase)

    data = {}

    for cId in configs:
        infile = scenario + '/' + phase + '/data/runtimes-' + cId + '.log'
        #accept ".log" or ".csv" files.
        if(not helper.isFile(infile)):
            infile = infile[:-3] + 'csv'
        if(not helper.isFile(infile)):
            print("Can't find the file " + infile + " or " + infile[:-3] + "log.")
            print("Exiting.")
            exit()
        data[cId] = {}        

        with open(infile) as f_in:
            if(verbose > 0):
                print("Loading data from: " + infile)

            for line in f_in:
                if('#' in line[0]):
                    continue
                items = line.split(',')
                inst = items[0].strip()
                seed = int(items[1].strip())
                size = int(items[2].strip())
                runtime = float(items[3])
                if(len(items) >= 4):
                    misc = items[4].strip()
                else:
                    misc = ''
                
                if(size not in data[cId].keys()):
                    data[cId][size] = {}
                if(inst not in data[cId][size].keys()):
                    data[cId][size][inst] = {}

                data[cId][size][inst][seed] = (runtime,misc)

                 
    helper.saveObj(scenario + '/' + phase,data,'data')

    return data
