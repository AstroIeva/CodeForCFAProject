import pandas as pd
import numpy as np
from matplotlib import pyplot as plt 
import matplotlib
import math
import statistics
import os

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
    
    #Obtaining the errors for counts
    countsErr = np.array(df.loc[:,"STAT_ERR"])
    
    return time, timeKs, counts, countsErr


def lightcurve(binsize,time,counts, countsErr):
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
    
    
    #Creating a list to store number of counts and err in each bin
    binCounts = []
    binCountsErr = []
    
    #Replace possible nan values to 0s to be able to sum with integers
    countsErr = np.nan_to_num(countsErr) 
    counts = np.nan_to_num(counts)
    
    countsInBin = 0 #Keep count for number of counts in each bin
    countsErrInBin = 0 #Keep count for error of counts in each bin
    indexCount = 0 #For knowing the corresponding time
    binIndex = 0
    prevTime = 0
    for i in counts:
        countsInBin = countsInBin + i
        countsErrInBin = countsErrInBin + countsErr[indexCount]
        if (time[indexCount]-prevTime)>binsize:
            binCounts.append(countsInBin)
            countsErrInBin = math.sqrt(countsErrInBin)
            binCountsErr.append(countsErrInBin)
            countsInBin = 0
            countsErrInBin = 0
            binIndex += 1
            prevTime = binsize*binIndex
        indexCount += 1 
    
    numBins = len(binCounts)
    #print(numBins)
    #print(binCounts)
    binTime = np.arange(time[0],time[-1],binsize)
    #print(len(binTime))
    #print(binTime)

    """
    #Obtain the number of bins:
    bins = time[-1]/binsize
    #Rounding bins up:
    numBins = math.ceil(bins)
    print(bins)
    print(numBins)

    #Creating a list to store number of counts and err in each bin
    binCounts = np.array([])
    binCountsErr = np.array([])

    countsInBin = 0 #Keep count for number of counts in each bin
    countsErrInBin = 0 #Keep count for error of counts in each bin
    for i in range(numBins):
        for j in range(len(time)):
            if i == 0 and time[j]<=binsize: #Counting the first bin
                print(i)
                print(time[j])
                countsInBin = countsInBin + counts[j]
                countsErrInBin = countsErrInBin + countsErr[j]**2
                print(countsErrInBin)
            else:
                if i != 0 and time[j]<(i+1)*binsize and time[j]>i*binsize: #Search for in between the binsize
                    countsInBin = countsInBin + counts[j]
                    countsErrInBin = countsErrInBin + countsErr[j]**2
        
        countsErrInBin = math.sqrt(countsErrInBin)
        #Update the arrays
        binCounts = np.append(binCounts, countsInBin) 
        binCountsErr = np.append(binCountsErr, countsErrInBin) 
        #Restart counting for next bin
        countsInBin = 0
        countsErrInBin = 0
        
    """     
            
    #Append with the last measurement twice to be able to draw the last bin line   
    #binCounts = np.append(binCounts, binCounts[-1])
    #binCountsErr = np.append(binCountsErr, binCountsErr[-1]) 

    
    
    return binTime[:-1],binCounts, numBins, binCountsErr


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

        
def plotting(time, timeKs, counts, countsErr, name):
    #List of binsizes in s:
    binsizes = [50, 100, 500, 1000, 1500]
      
    #Creating multiple plots
    fig, axes = plt.subplots(5, 2, figsize=(20, 20))
    
    #Naming the plots
    fig.suptitle("Plots for " + str(name[:-12]), fontsize=25)
    #Setting the plot size
    matplotlib.rcParams.update({'font.size': 20})
    
    for i in range(5):
        #Make calculations for lightcurve plot
        binTime,binCounts, numBins, binCountsErr = lightcurve(binsizes[i],time,counts, countsErr)
        
        #Plot for lightcurve
        binsize = binsizes[i]
        print(binsize)
        print(numBins)
        axes[i,0].set_title(f"Lightcurve with {numBins} bins, and binsize = {binsize}s")
        axes[i,0].set_xlabel("Time (ks)") 
        axes[i,0].set_ylabel("Total counts") 
        axes[i,0].step([0,binsizes[i]/2/1000],[binCounts[0],binCounts[0]], where="pre", color = "blue") #Draw the first bit to the plot
        axes[i,0].step([((binTime[-1]+binsizes[i]/2)/1000),((binTime[-1]+binsizes[i])/1000)],[binCounts[-1],binCounts[-1]], where="post", color = "blue") #Draw the last bit to the plot
        axes[i,0].step((binTime+binsizes[i]/2)/1000,binCounts, where="mid", color = "blue")
        axes[i,0].errorbar((binTime+binsizes[i]/2)/1000, binCounts, yerr=binCountsErr, fmt = "o" ,ecolor="black",linewidth = .5, 
                     capsize = 1,color="black")
        
        
        #Make calculations for running average plot
        timeAvg,avgBinCount,avgBinCountErr = runningAvg(binTime, binCounts)
        
        #Plot for counts in 3 bins avg
        axes[i,1].set_title(f"Running average over 3 bins for binsize = {binsize}s")
        axes[i,1].set_xlabel("Time (ks)") 
        axes[i,1].set_ylabel("Counts avg") 
        axes[i,1].plot(timeAvg/1000,avgBinCount, color = "crimson", linewidth=2.5) 
        axes[i,1].errorbar(timeAvg/1000,avgBinCount,yerr = avgBinCountErr, fmt = '.', color = "black",linewidth = .5,
                 capsize = 1) 
    
    fig.tight_layout()
    plt.savefig(name[:-12]+".png")
    #plt.show()"""
        
def main():
    
    #Folder name (must be same path where the code is) 
    dirPath = os.getcwd()
    
    # list to store files
    fileNames = []
    #Iterate directory
    for path in os.listdir(dirPath):
        # check if current path is a file
        if os.path.isfile(os.path.join(dirPath, path)):
            #Look for txt files only:
            if path[-3:]=="txt":
                if path != "README.txt":
                    fileNames.append(path)
            
    print(fileNames)
    for name in fileNames: #Run through all data files in directory
        print(name)
        #Read data
        time, timeKs, counts, countsErr = getData(name)
        
        #Process it and plot
        plotting(time, timeKs, counts, countsErr, name)
        
main()



