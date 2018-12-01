#!/usr/bin/env python
from scipy.stats import norm
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt

x1 = np.loadtxt('../../results/phase_1/error_2000.dat.0', usecols=(4,))
x2 = np.loadtxt('../../results/phase_1/error_100.dat.0', usecols=(4,))

# Empirical average and variance are computed
avg1 = np.mean(x1)
var1 = np.var(x1)
avg2 = np.mean(x2)
var2 = np.var(x2)
# From that, we know the shape of the fitted Gaussian.
pdf_x1 = np.linspace(np.min(x1),np.max(x1),100)
pdf_y1 = 1.0/np.sqrt(2*np.pi*var1)*np.exp(-0.5*(pdf_x1-avg1)**2/var1)

pdf_x2 = np.linspace(np.min(x2),np.max(x2),100)
pdf_y2 = 1.0/np.sqrt(2*np.pi*var2)*np.exp(-0.5*(pdf_x2-avg2)**2/var2)

# Then we plot :
plt.figure()
#plt.hist(data,30,normed=True)
plt.plot(pdf_x1,pdf_y1,'k--')
plt.plot(pdf_x2,pdf_y2,'k--')
#plt.legend(("Fit","Data"),"best")
#plt.show()

#(mu, sigma) = norm.fit(x1)
# the histogram of the data
n, bins, patches = plt.hist(x1, 100, normed=True, facecolor='green')
n2, bins2, patches2 = plt.hist(x2, 100, normed=True, facecolor='red')

# add a 'best fit' line
#y = mlab.normpdf( bins, mu, sigma)
#l = plt.plot(bins, y, 'r--', linewidth=1)

plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
#plt.axis([40, 160, 0, 0.03])
plt.grid(True)

plt.show()
