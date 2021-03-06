#!/usr/bin/env python
import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#mu, sigma = 100, 15
filename='/export/home/klast308/wc/results/25/phase_1/error_2000.dat.1'
data = np.loadtxt(filename,usecols=(1,2,3),skiprows=1)
idstr = os.path.join(os.path.basename(os.path.dirname(filename)),os.path.basename(filename))

pred=data[:,1]
real=data[:,2]

plt.rc('xtick', labelsize=18)
plt.rc('ytick', labelsize=18) 
plt.rcParams.update({'font.size': 20})
plt.rcParams['axes.linewidth'] = 2

plt.figure(1)
plt.plot(pred,real,'.')
plt.xlabel('Predicted [eV]')
plt.ylabel('Real [eV]')
plt.title('')
plt.grid(True)

plt.show()
