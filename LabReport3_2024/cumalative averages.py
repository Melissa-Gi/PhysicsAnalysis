import math
import numpy as np 
import scipy as sp
import matplotlib.ticker as ticker
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

channelData = np.zeros(1024)
countsData = np.zeros(1024)
udata = np.zeros(1024)
i = 0

def readData(filename):
        f = open(filename,'r')         # open file
        for y in range (0,22):
                f.readline()                   # read and ignore header
        global i
        i = 0
        global channelData
        global countsData
        global udata

        for line in f:                          # loop over lines
                line.strip()
                columns = line.split(',')
                channelData[i] = float(columns[0])
                countsData[i] = float(columns[1])
                i = i + 1
        f.close()

filenames = ['MOB/30CPI_1CM_20MS.csv','MOB/10CPI_4_4CM_20MS.csv','MOB/1CPI_22CM_20MS.csv','MOB/4CPI_11_4CM_20MS.csv','MOB/30CPI_1CM_40MS.csv','MOB/30CPI_1CM_100MS.csv','MOB/30CPI_1CM_200MS.csv']
labels = ['20 ms, 30 avg', '20ms, 10 avg','20ms, 1 avg','20ms, 4 avg','40 ms','100ms','200ms']
for i in range(0,len(labels)):
        fig, ax = plt.subplots()

        ax.set_xlabel("Sequence number j of time intervals of the count series")
        ax.set_ylabel("Cumulative average, r(j)")
        ax.set_xlim(0,1050)
        
        readData(filenames[i])        
        
        dwell_time = 0.02
        intervals=np.arange(1,1025,1)

        cumulative_average = np.zeros(1024)

        std_dev = np.std(countsData)
        for i in range(0,1024):
                cumulative_average[i] = sum(countsData[0:i+1])/(i+1)
                udata[i]=std_dev/np.sqrt(channelData[i]+1)

        msg = 'Uncertainty related to the number of counts, N: {}/{}'.format(r'$\sigma$', r'$\sqrt{N}$')
        ax.fill_between(intervals, cumulative_average-udata, cumulative_average+udata, alpha=.5, linewidth=0,label=msg)
        plt.plot(intervals,cumulative_average,label='Cumulative average of counts')
        msg = r'$\mu$=' + str(round(cumulative_average[i],3))

        plt.scatter(intervals[i],cumulative_average[i],label=msg)
        #ax.set_ylim(0,np.max(cumulative_average+udata)+0.2*np.max(cumulative_average+udata))
        ax.legend()
        plt.show()

