import math
import statistics
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

def poisson(mu,r):
        return (np.power(mu,r))*np.exp(-mu)/sp.special.factorial(r)

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
mus = [32.10033333333333,8.879733333333334,0.8551,2.595266666666667,64.20066666666666,160.50166666666667,321.00333333333333]
for x in range(0,7):
        readData(filenames[x])

        fig, ax = plt.subplots()

        ax.set_xlabel("Number of counts per counting interval")
        ax.set_ylabel("Frequency F")
        bins = np.arange(np.min(countsData),np.max(countsData)+1)
        mean = np.mean(countsData)
        
        varianceHist = np.var(countsData)
        msg = labels[x] + ': '+ r'$\sigma^2$'+ ' of binned data =' + str(round(varianceHist,2))
        plt.hist(countsData, bins=bins-0.5,histtype='bar',color='lightsteelblue',ec='black',label=msg)
        plt.axvline(mean,color='blue',label=('Sample mean ' + r'$\mu$' + '= '+str(round(mean,2))))

        #Uncertainties:
        frequency, binEdges=np.histogram(countsData,bins=bins)
        msg='Uncertainty '+ r'$\sqrt{F}$'
        uncert=np.round(np.sqrt(frequency),decimals=2)
        plt.errorbar(bins[0:-1], frequency, yerr=uncert, label=msg, capsize=2,linewidth=0.7,linestyle='none')
        
        poisDist = stats.poisson.pmf(bins,mus[x])
        poisDist = poisDist*1048

        dataset=[]
        for y in range(len(poisDist)):
                for z in range(0,round(poisDist[y])):
                        dataset.append(bins[y])
        variancePois = statistics.variance(dataset)             
        
        msg = 'Poisson distribution with '+ r'$\sigma^2$'+ '= ' + str(round(variancePois,2))

        plt.plot(bins,poisDist,color='red',label=msg)

        ax.set_ylim(0,np.max(frequency)+10)
        
        if x > 3 or x==0:
                sd = np.std(countsData)
                gaussian = stats.norm.pdf(bins,mus[x],sd)*1048
                plt.plot(bins,gaussian,color='green',label='Gaussian using population mean') 

                difference = np.abs(poisDist[0:-1]-frequency)
                percent = 100*np.sum(difference)/np.sum(poisDist[0:-1])
                print('Largest percent deviation for Poisson',str(percent),'%')
                
                difference = np.abs(gaussian[0:-1]-frequency)
                percent = 100*np.sum(difference)/np.sum(gaussian[0:-1])
                print('Largest percent deviation for Gaussian',str(percent),'%')
                
                print()
                
        plt.legend()
        plt.show()


