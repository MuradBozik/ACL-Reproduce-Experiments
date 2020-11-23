#Author: Yasha Pushak
#Created: Some time in 2015
#A collection of general-purpose helper functions I have found useful over the years.

import time
import math
import pickle
import string, random
from contextlib import contextmanager
import os
import numpy as np
import glob
import datetime

def generateID(size=6, chars=string.ascii_uppercase + string.digits):
    #generate a random ID for identifying SMAC runs
    return ''.join(random.choice(chars) for _ in range(size))


def randSeed():
    return random.randint(0,2147483647)


def mkdir(dir):
    #Author: Yasha Pushak
    #Last updated: January 3rd, 2017
    #An alias for makeDir
    makeDir(dir)

def makeDir(dir):
    #Only creates the specified directory if it does not already exist.
    #At some point it may be worth adding a new feature to this that saves the old directory and makes a new one if one exists.
    if(not os.path.isdir(dir)):
        os.system('mkdir '  + dir)


def isDir(dir):
    #Author: Yasha Pushak
    #last updated: March 21st, 2017
    #Checks if the specified directory exists.
    return os.path.isdir(dir)


def isFile(filename):
    #Author: Yasha Pushak
    #Last updated: March 21st, 2017
    #CHecks if the specified filename is a file.
    return os.path.isfile(filename)


#Code taken from http://stackoverflow.com/questions/19201290/how-to-save-a-dictionary-to-a-file

def saveObj(dir, obj, name ):
    with open(dir + '/'+ name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def loadObj(dir, name ):
    with open(dir + '/' + name + '.pkl', 'rb') as f:
        return pickle.load(f)


def calStatistic( lst, statistic, cutoff=float('inf') ):
    if statistic.lower() == 'par10':
        statistic = "mean"
        lst = np.array(lst)
        lst[np.where(lst >= cutoff)] = cutoff*10
        lst = list(lst)
    if statistic.lower() == "mean":
        return sum( lst )/len(lst)
    if statistic.lower() == "median":
        statistic = "q50"
    if statistic[0].lower() == "q":
        percent = float( statistic[1:] )
        if percent<1:
            percent *= 100

        lst = sorted(lst)
        I = len(lst)*(percent/100.0) - 1
        #Check if I is an integer
        if(int(I) - I == 0):
           return (lst[int(I)] + lst[int(I+1)])/2
        else:
           return lst[int(math.ceil(I))]

        #YP: The original code here used to always return a lower-bound on
        #the quantiles. I have changed this..
        #This is what Zongxu used to have: 
        #return sorted(list)[ int(len(list)*percent/100)-1 ]
    raise ValueError('Invalid summary statistic input: ' + statistic)



@contextmanager
def acquireLock(lockFile,delay=10):
    #Author: YP
    #Last updated: 2018-08-21
    #A crudge, text-file based lock system.
    #Note guaranteed to always work. Out of several
    #hundred runs I saw it fail once.
    #If things get wierd, delete the "lock file" and
    #kill the processes, then try again.
    myId = getLock(lockFile,delay)
    try:
        yield
    finally:
        releaseLock(lockFile,myId)



def getLock(lockFile,delay=10):
    locked = True
    myId = generateID(100)
    while(locked):
      try:
        lockState = 'unlocked'
        if(isFile(lockFile)):
            #Check the state of the lock
            with open(lockFile) as f_in:
                lockState = f_in.read()
        if(lockState == 'unlocked'):
            #The lock is unlocked, try to acquire it
            with open(lockFile,'w') as f_out:
                f_out.write(myId)
            #We might have attempted to acquire the lock at the same time as someone else causing a race condition. 
            #Go to sleep for a bit and then check if the lock state was successfully updated
            time.sleep(random.randrange(delay-1,delay+1))
            with open(lockFile) as f_in:
                lockState = f_in.read()
            if(lockState == myId):
                #Lock acquired!
                return myId
            elif(not len(lockState) == 100 and not lockState == 'unlocked'):
                #It seems like someone else tried to get the lock at the same time and a race condition occured.
                #We will reset the state to unlocked, then go to sleep. (So the other person who tried to acquire the lock might get it this time)
                with open(lockFile,'w') as f_out:
                    f_out.write('unlocked')
      except: #In case a race condition causes an exception.
          print("An exception was caught waiting for the lock.")
      print("Waiting for lock " + lockFile + "...")
      time.sleep(random.randrange(delay-1,delay+1))


def releaseLock(lockFile,myId):
    if(isFile(lockFile)):
        with open(lockFile) as f_in:
            lockState = f_in.read()
        if(not lockState == myId):
            print("WARNING: I just tried to release a lock but it was not mine to release...")
        else:
            with open(lockFile,'w') as f_out:
                f_out.write('unlocked')
    else:
        print("WARNING: I was told to delete a lock that does not exist...")



