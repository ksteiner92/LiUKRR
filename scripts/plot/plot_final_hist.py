#!/usr/bin/env python
import os
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

#mu, sigma = 100, 15
filename='/export/home/klast308/wc/results/phase_1/error_3000.dat.2'
data = np.loadtxt(filename,usecols=(1,2,3),skiprows=1)
idstr = os.path.join(os.path.basename(os.path.dirname(filename)),os.path.basename(filename))

pred=data[:,1]
real=data[:,2]

plt.figure(1)
n, bins, patches = plt.hist(pred-real, 500, normed=1, facecolor='green', alpha=0.75)
plt.xlabel('Error')
plt.ylabel('Distribution')
plt.title(idstr)
plt.grid(True)
plt.xlim([-100,100])

plt.figure(2)
n, bins, patches = plt.hist(pred, 100, normed=1, facecolor='green', alpha=0.75)
plt.xlabel('Prediction')
plt.ylabel('Distribution')
plt.title(idstr)
plt.grid(True)

plt.show()
