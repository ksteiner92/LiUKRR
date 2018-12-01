#!/usr/bin/env python
import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#mu, sigma = 100, 15
filename='/export/home/klast308/wc/results/25/phase_all/error_3000.dat.2'
data = np.loadtxt(filename,usecols=(1,2,3),skiprows=1)
idstr = os.path.join(os.path.basename(os.path.dirname(filename)),os.path.basename(filename))

pred=data[:,1]
real=data[:,2]

plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18) 
plt.rcParams.update({'font.size': 20})
plt.rcParams['axes.linewidth'] = 2

plt.figure(1)
n, bins, patches = plt.hist(pred-real, 150, normed=1, facecolor='green', alpha=0.75)
plt.xlabel('Error [eV]')
plt.ylabel('Distribution')
plt.title('')
plt.grid(True)
plt.xlim([-100,100])

plt.figure(2)
n, bins, patches = plt.hist(pred, 75, normed=1, facecolor='green', alpha=0.75, label='predicted data')
plt.hist(real, 75, normed=1, facecolor='red', alpha=0.5, label='real data')
plt.xlabel('Values [eV]')
plt.ylabel('Distribution')
plt.title('')
plt.grid(True)
plt.legend(loc=2,fontsize=16)
plt.show()
