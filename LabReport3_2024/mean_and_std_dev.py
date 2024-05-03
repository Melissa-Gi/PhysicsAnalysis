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

def f(x, A, mu, sig):
    return A*np.exp(-(x-mu)**2/(2*sig**2))

def readData(filename):
        f = open(filename,'r')         # open file
        for y in range (0,22):
                f.readline()                   # read and ignore header
        global i
        i = 0
        global channelData
        global countsData

        for line in f:                          # loop over lines
                line.strip()
                columns = line.split(',')
                channelData[i] = float(columns[0])
                countsData[i] = float(columns[1])
                udata[i] = 1
                i = i + 1
        f.close()

filenamesA = ['MOB/4_4CM.csv','MOB/11_4CM.csv','MOB/22CM.csv','MOB/1CM.csv']
namesA = ['10 counts long spec','4 counts long spec','1 counts long spec','30 counts long spec']
for x in range(0, 4):
    readData(filenamesA[x])
    print(namesA[x])
    print('Mean:',np.mean(np.sum(countsData)/30000))
    print('Standard deviation:',np.std(countsData/30000))
    print()

namesA = ['30 counts long spec scaled for 40ms','30 counts long spec scaled for 100ms','30 counts long spec scaled for 200ms']
scales=[15000,6000,3000]
for x in range(0,3):
    print(namesA[x])
    print('Mean:',np.mean(np.sum(countsData)/scales[x]))
    print('Standard deviation:',np.std(countsData[205:310]/scales[x]))
    print()


filenames = ['MOB/10CPI_4_4CM_20MS.csv','MOB/4CPI_11_4CM_20MS.csv','MOB/1CPI_22CM_20MS.csv','MOB/30CPI_1CM_20MS.csv']
names = ['Avg 10 counts','Avg 4 counts','Avg 1 count','Avg 30 counts']
for x in range(0, 4):
    print()
    readData(filenames[x])
    print(names[x])
    print('The mean is:',np.mean(countsData),'\nThe standard deviation is:',np.std(countsData))

    
filenames = ['MOB/30CPI_1CM_40MS.csv','MOB/30CPI_1CM_100MS.csv','MOB/30CPI_1CM_200MS.csv']
names = ['40 ms','100 ms','200 ms']
for x in range(0,3):
    readData(filenames[x])
    print(names[x])
    print('The mean is:',np.mean(countsData),'\nThe standard deviation is:',np.std(countsData))
    
