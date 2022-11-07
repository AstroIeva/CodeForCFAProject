import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
import math
import statistics
import os

"""
==================================================
                    Edit here
==================================================
"""
#Folder name 
dirPath = "/home/ievajankute/Documents/Work/Project"   

#binsize
binsize = 1000

"""
==================================================
"""

def getData(name):
    
    """
    This function reads from a txt file the time and counts.
    
    Parameters
    ----------
    
    name: name of the txt file to read
    
    Returns
    -------
    
    time: time list of the measurements in seconds starting at 0
    
    timeKs: time list in kiloseconds
    
    counts: counts list
    
    """
    #Reading the file
    df = pd.read_csv(name,sep=" ")

    """This file has the following headers:
        TIME_BIN
        TIME_MIN
        TIME
        TIME_MAX
        COUNTS
        STAT_ERR
        AREA
        EXPOSURE
        COUNT_RATE
        COUNT_RATE_ERR
        """
    
    #Obtaining the time measurement
    timeRead = np.array(df.loc[:,"TIME"])

    #Setting the time to be from the start of the measurement
    time = timeRead-timeRead[0]

    #Convert time from s to ks
    timeKs = time/1000

    #Obtaining the counts measurement
    counts = np.array(df.loc[:,"COUNTS"])
    return time, timeKs, counts


def cumulativeCountsCalc(counts):
    """
    This function calculates the cumulative counts.
    

    Parameters
    ----------
    counts: Counts list

    Returns
    -------
    cumulativeCounts: A list of count increase over time

    """
    #Obtaining cumulative counts data (increase in counts over time)
    cumulativeCounts = np.array([]) #New array to store the data
    totCounts = 0 #Total counts up to this point in time 
    for i in range(len(counts)):
        totCounts = totCounts + counts[i] #Increase in total count
        cumulativeCounts = np.append(cumulativeCounts, totCounts) #Update the array
        
    return cumulativeCounts


def lightcurve(binsize,time,counts):
    """
    This function makes bins of the chosen binsize and adds counts to the corresponding time size bin.

    Parameters
    ----------
    binsize : Size of the bin in seconds
    time : time list as read from data
    counts : counts list as read from data

    Returns
    -------
    binTime : time list for plotting the lightcurve
    binCounts : binned counts list for the lightcurve
    numBins : number of bins used

    """
    """Here we also round up the bin size, so the last bin has contribution from less time"""
    
    #Lightcurve data
    #Obtain the number of bins:
    bins = time[-1]/binsize
    #Rounding bins up:
    numBins = math.ceil(bins)

    #Creating a list to store number of counts in each bin
    binCounts = np.array([])

    countsInBin = 0 #Keep count for number of counts in each bin
    for i in range(numBins):
        for j in range(len(time)):
            if i == 0 and time[j]<binsize: #Counting the first bin
                countsInBin = countsInBin + counts[j]
            else:
                if time[j]<(i+1)*binsize and time[j]>i*binsize: #Search for in between the binsize
                    countsInBin = countsInBin + counts[j]
        binCounts = np.append(binCounts, countsInBin) #Update the array
        countsInBin = 0 #Restart counting for next bin
            
    #Append with the last measurement twice to be able to draw the last bin line   
    binCounts = np.append(binCounts, binCounts[-1])

    binTime = np.arange(time[0],time[-1]+binsize,binsize)
    
    return binTime,binCounts, numBins


def runningAvg(binTime, binCounts):
    
    """
    This function takes the number of counts in middle back and front, and averages that with also calculating 
    stdev for errors.
    
    Parameters
    ----------
    binTime : binned time list as for lightcurve plot
    binCounts : binned number of counts as for lightcurve plot
    
    Returns
    -------
    timeAvg : time list, where time indicated corresponds to avg counts at the middle
    avgBinCount : average bin counts list
    avgBinCountErr : average bin count list error using stdev.
    """
    
    #Making counter around which to count bins:
    binCountCenter = 1
    avgOver = 3 #Average over what number
    #Creating arrays to store averages and their errors:
    avgBinCount = np.array([])
    avgBinCountErr = np.array([])
    timeAvg = np.array([])
    summed = 0 #Counter for sum
    i = 0
    current3avgs = np.array([]) #List of current 3 bins for finsing error

    while i < len(binCounts)-(avgOver-1):#Run through all values in bin counts
        for j in range(avgOver):
            if i<=binCountCenter+1:
                summed = summed + binCounts[i]
                current3avgs = np.append(current3avgs, binCounts[i])
                i = i + 1
        #Finding time:
        timeAvg = np.append(timeAvg, binTime[binCountCenter])
        #Finding error:
        err = statistics.stdev(current3avgs)
        #Finding average:
        avg = summed/avgOver
        #Updating lists:
        avgBinCount = np.append(avgBinCount, avg)
        avgBinCountErr = np.append(avgBinCountErr, err)
        #Reset and continue with updated loop values
        current3avgs = np.array([])
        summed = 0 
        binCountCenter = binCountCenter + 1
        i = binCountCenter  - 1
    
    return timeAvg,avgBinCount,avgBinCountErr
    

def plotting(timeKs,cumulativeCounts, counts, binTime,binCounts, avgBinCount, avgBinCountErr,timeAvg,numBins,binsize,name):
    """
    Function for creating and plotting all the neccessary plots.

    Parameters
    ----------
    timeKs : time list in kiloseconds
    cumulativeCounts : A list of count increase over time
    counts : counts list as read from txt file
    binTime : binned time list as for lightcurve plot
    binCounts : binned number of counts as for lightcurve plot
    avgBinCount : average bin counts list
    avgBinCountErr : average bin count list error using stdev.
    timeAvg : time list, where time indicated corresponds to avg counts at the middle
    numBins : number of bins used for lightcurve
    binsize : Size of the bin in seconds
    name : name of the txt file

    Returns
    -------
    None.

    """

    #Creating multiple plots
    fig, axes = plt.subplots(2, 2, figsize=(15, 15))

    fig.suptitle("Plots for " + str(name[:-12]), fontsize=18)
    
    #Making cumulative counts plot
    plt.subplot(2, 2, 1) #Plot in the top left for lightcurve
    plt.title(f"Lightcurve with {numBins} bins, and binsize = {binsize}s")
    plt.xlabel("Time (ks)") 
    plt.ylabel("Total counts per bin") 
    plt.step(binTime/1000,binCounts, where="post")

    plt.subplot(2, 2, 2) #Plot in the top right for counts in 3 bins avg
    plt.title("Running average over 3 bins")
    plt.xlabel("Time (ks)") 
    plt.ylabel("Total counts average per 3 bins") 
    plt.plot(timeAvg/1000,avgBinCount, color = "crimson", linewidth=2.5) 
    plt.errorbar(timeAvg/1000,avgBinCount,yerr = avgBinCountErr, fmt = '.', color = "black",linewidth = .5,
             capsize = 1) 

    plt.subplot(2, 1, 2) #Plot in the bottom for cumulative counts
    plt.title("Cumulative counts vs time")
    plt.xlabel("Time (ks)") 
    plt.ylabel("Total counts") 
    plt.plot(timeKs,cumulativeCounts, color = "purple")

    fig.tight_layout()
    plt.savefig(name[:-12]+".png")
    plt.show()

    

def main(dirPath,binsize):
    
    # list to store files
    txtFileNames = []

    # Iterate directory
    for path in os.listdir(dirPath):
        # check if current path is a file
        if os.path.isfile(os.path.join(dirPath, path)):
            #Look for txt files only:
            if path[-3:]=="txt":
                txtFileNames.append(path)
    
    for name in txtFileNames: #Run through all data files in directory
        #Read data
        time, timeKs, counts = getData(name)
    
        #Make calculations for cumulative counts plot
        cumulativeCounts = cumulativeCountsCalc(counts)
    
        #Make calculations for lightcurve plot
        binTime,binCounts, numBins = lightcurve(binsize,time,counts)
    
        #Make calculations for running average plot
        timeAvg,avgBinCount,avgBinCountErr = runningAvg(binTime, binCounts)
    
        #Making all the plots
        plotting(timeKs,cumulativeCounts,counts, binTime,binCounts, avgBinCount, avgBinCountErr, timeAvg,numBins,binsize, name)

            
main(dirPath, binsize)




