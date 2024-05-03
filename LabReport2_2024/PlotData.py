import math
import numpy as np 
import scipy as sp
import matplotlib.ticker as ticker
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d

channelData = np.zeros(2048)
countsData = np.zeros(2048)
udata = np.zeros(2048)
i = 0

def f(x, A, mu, sig):
    return A*np.exp(-(x-mu)**2/(2*sig**2))
def decay(t,N0,lamda,B):
    return N0*np.exp(-lamda*t) + B
def convertCtoE(channels):
    return (channels+6.858128897802162)/0.26847222953379374

def readData(filename):
        f = open(filename,'r')         # open file
        for y in range (0,22):
                f.readline()                   # read and ignore header
        global i
        i = 0
        global channelData
        global countsData

        for line in f:                          # loop over lines
                columns = line.split('\t')
                channelData[i] = float(columns[0])
                countsData[i] = float(columns[2])
                udata[i] = 1
                i = i + 1
        f.close()

#----------------------Plot to show gaussian fits and raw data----------------------#
#plotting variables
channels = convertCtoE(np.arange(0,2048,1))
bins = convertCtoE(np.arange(0,2049,1))
time = np.linspace(0,3600,120)
totalSpectra = np.zeros(2048)

decayCount_bottom = np.zeros(120)
decayCount_middle = np.zeros(120)
decayCountAl = np.zeros(120)
udataAl = np.zeros(120)
decayCountMg = np.zeros(120)
udataMg = np.zeros(120)

half_life_of_Al = 2.245*60
lambaAl = np.log(2)/half_life_of_Al
lambaSi = 1000000

fig1, ax1 = plt.subplots()

ax1.set_xlabel("Energy (keV)")
ax1.set_ylabel("Count #")
ax1.xaxis.set_major_locator(ticker.MultipleLocator(200))
ax1.set_xlim(74,2900)

#---------Bottom----------#
#Sum up spectra
for j in range(1,121):
    filename = 'trial3_bottom/trial3_' + str(j) + '.tsv'
    readData(filename)
    totalSpectra = totalSpectra + countsData
    
    decayCount_bottom[j-1] = sum(countsData)
    decayCountAl[j-1] = sum(countsData[461:492])
    udataAl[j-1] = np.sqrt(sum(countsData[450:500]))
    decayCountMg[j-1] = sum(countsData[212:231])
    udataMg[j-1] = np.sqrt(sum(countsData[209:231]))

both_trial_total = totalSpectra
bottom_total = totalSpectra
channelData = convertCtoE(channelData)
udata = 0.0293*channelData + 8.3983

plt.stairs(bottom_total,bins,label='Gamma Spectrum of The Bottom Cylinder')

#Zero variable containers
totalSpectra = np.zeros(2048)
channelData = np.zeros(2048)
countsData = np.zeros(2048)

#---------middle----------#
#Sum up spectra
for j in range(1,121):
    filename = 'trial1_middle_cylinder/trial1_' + str(j) + '.tsv'
    readData(filename)
    totalSpectra = totalSpectra + countsData
    
    decayCount_middle[j-1] = sum(countsData)

both_trial_total = both_trial_total + totalSpectra
middle_total = totalSpectra
channelData = convertCtoE(channelData)

#plt.stairs(totalSpectra,bins,color='dimgrey',label='Middle Cylinder')

#------------Fit and process of data-----------#

#Peak 1
A = 1985
mu = 853
sig = 80
p0 = [A, mu, sig]   #List of initial guesses
channelModel2 = np.linspace(792,912,1000)
popt1,pcov1 = curve_fit(f, channels[205:238], bottom_total[205:238],p0,sigma=udata[205:238],absolute_sigma = True)

print('\nPeak 1: Magnesium 27 Decay, 1st Excited Al27')
perr = np.sqrt(np.diag(pcov1))
print('Amplitude:',popt1[0],'+/-',perr[0],'\nCentroid:',popt1[1],'+/-',perr[1],'keV','\nStandard deviation:',np.absolute(popt1[2]),'+/-',perr[2],'keV')
y2fit=f(channelModel2, *popt1) #Vectorised
plt.plot(channelModel2, y2fit, color='red',linewidth = 1)

#Peak 2
A = 301
mu = 1368.626
sig = 100
p0 = [A, mu, sig]   #List of initial guesses
channelModel2 = np.linspace(1355,1402,1000)
popt1,pcov1 = curve_fit(f, channels[357:369], bottom_total[357:369],p0,sigma=udata[357:369],absolute_sigma = True)

print('\nPeak 2: Isotope of Sodium: Possibly Na 24 producing Mg 24 through Beta- decay')
perr = np.sqrt(np.diag(pcov1))
print('Amplitude:',popt1[0],'+/-',perr[0],'\nCentroid:',popt1[1],'+/-',perr[1],'keV','\nStandard deviation:',popt1[2],'+/-',perr[2],'keV')
y2fit=f(channelModel2, *popt1) #Vectorised
plt.plot(channelModel2, y2fit, color='red',linewidth = 1)

#Peak 3
A = 326
mu = 1475
sig = 100
p0 = [A, mu, sig]   #List of initial guesses
channelModel3 = np.linspace(1426,1530,1000)
popt1,pcov1 = curve_fit(f, channels[376:404], bottom_total[376:404],p0,sigma=udata[376:404],absolute_sigma = True)

print('\nPeak 3: Isotope of Sodium: Possibly Na 28 producing Mg 28 through Beta- decay')
perr = np.sqrt(np.diag(pcov1))
print('Amplitude:',popt1[0],'+/-',perr[0],'\nCentroid:',popt1[1],'+/-',perr[1],'keV','\nStandard deviation:',popt1[2],'+/-',perr[2],'keV')
y3fit=f(channelModel3, *popt1) #Vectorised
plt.plot(channelModel3, y3fit, color='red',linewidth = 1)

#Peak 4
A = 359
mu = 1800
sig = 100
p0 = [A, mu, sig]   #List of initial guesses
channelModel4 = np.linspace(1705,1893,1000)
popt1,pcov1 = curve_fit(f, channels[450:500], bottom_total[450:500],p0,sigma=udata[450:500],absolute_sigma = True)

print('\nPeak 3: Aluminium 28 Decay')
perr = np.sqrt(np.diag(pcov1))
print('Amplitude:',popt1[0],'+/-',perr[0],'\nCentroid:',popt1[1]-10,'+/-',perr[1],'keV','\nStandard deviation:',popt1[2],'+/-',perr[2],'keV')
y4fit=f(channelModel4, *popt1) #Vectorised
plt.plot(channelModel4, y4fit, color='red', label = 'Gaussian fits for unique peaks of interest',linewidth = 1)

plt.yscale('log')
plt.legend()
plt.show()

#----------------------Decay Curve----------------------#
print('\nDecay Curve')

fig, ax = plt.subplots()
ax.set_xlabel("Time (s)")
ax.set_ylabel("Count Rate (counts per 30 second interval)")
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.2))

ax.set_xlim(0,2500)
#decayCountAl_a = decayCountAl_a/30
#decayCountAl_b = decayCountAl_b/30

#Plot of all decays for middle and bottom isotopes
bins = np.arange(0,3630,30)

plt.stairs(decayCountAl,bins, label='Al28 data from bottom cylinder')
ax.errorbar(bins[0:120]+15,decayCountAl,udataAl,capsize=1,color='tab:blue',linestyle='None',linewidth=0.5)

ax.errorbar(bins[0:120]+15,decayCountMg,udataMg,capsize=1,color='tab:orange',linestyle='None',linewidth=0.5)
plt.stairs(decayCountMg,bins, label='Mg27 data from bottom cylinder')

plt.stairs(decayCount_bottom,bins, label='Total decay of bottom cylinder')
#plt.stairs(decayCount_middle,bins, color='dimgrey',label='Total decay of middle cylinder')

#Aluminium decay fit
background = np.average(decayCountAl[43:119])
p0 = [decayCountAl[0],lambaAl,background]   #List of initial guesses
popt1,pcov1 = curve_fit(decay, time, decayCountAl,p0,sigma=udataAl,absolute_sigma = True)
print('\nDecay Curve of Aluminium')
perr = np.sqrt(np.diag(pcov1))
print('Half life =', np.log(2)/popt1[1]/60,'minutes')
print('N0:',popt1[0],'+/-',perr[0],'\nLambda:',popt1[1],'+/-',perr[1],'\nBackground energy:',popt1[2],'+/-',perr[2])
y1fit=decay(time, *popt1) #Vectorised
plt.plot(time, y1fit, '-b', label = 'Exponential fit for Al28 Decay',linewidth = 1)

#Magnesium decay fit
background = np.average(decayCountMg[43:119])
p0 = [decayCountMg[0],lambaAl,background]   #List of initial guesses
popt1,pcov1 = curve_fit(decay, time, decayCountMg,p0,sigma=udataMg,absolute_sigma = True)
print('\nDecay Curve of Magnesium')
perr = np.sqrt(np.diag(pcov1))
print('Half life =', np.log(2)/popt1[1]/60,'minutes')
print('N0:',popt1[0],'+/-',perr[0],'\nLambda:',popt1[1],'+/-',perr[1],'\nBackground energy:',popt1[2],'+/-',perr[2])
y1fit=decay(time, *popt1) #Vectorised
plt.plot(time, y1fit, '-r', label = 'Exponential fit for Mg27 Decay',linewidth = 1)

plt.yscale('log')
plt.legend()
plt.show()
