#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

IQR = lambda x        : st.scoreatpercentile(x, 75.0) - st.scoreatpercentile(x, 25.0)

def nbins(sample, range_ = None) :
	if range_ is None : mn, mx = sample.min(), sample.max()
	else              : mn, mx = range_

	mask    = (sample >= mn) & (sample <= mx)
	binsize = (2 * IQR(sample[mask]) / mask.sum() ** (1. / 3))

	return (mx - mn) / binsize, mn, mx, binsize

mu = np.genfromtxt("../inputs/run07_sfh_catalog.txt", missing = '""', usecols = [14])
Z  = np.genfromtxt("../inputs/run07_sfh_catalog.txt", missing = '""', usecols = [12])
Z1 = Z[Z<0.2]
Z2 = Z[Z>0.2]

plt.figure()
n, xi, xf, bs = nbins(Z1, range_ = (0.058920444050700002, 0.2))
counts, xbins = np.histogram(Z1, n, range = (xi, xf))
ibn = np.digitize(Z1, xbins)
idx = [np.where(ibn == i)[0] for i in xrange(1, xbins.size)]
Z1_idx = np.concatenate([bin[:counts.min()] for bin in idx])
plt.hist(Z1, n, histtype = "stepfilled", color = "#BFBFBF", lw = 0, range = (xi, xf), normed = True)
plt.hist(Z1[Z1_idx], n, histtype = "step", color = "#FFA500", lw = 2, range = (xi, xf), normed = True)
plt.plot([0.06, 0.2], [1/0.14,1/0.14], "-r", lw = 2)

n, xi, xf, bs = nbins(Z2, range_ = (0.2, 2.5))
counts, xbins = np.histogram(Z2, n, range = (xi, xf))
ibn = np.digitize(Z2, xbins)
idx = [np.where(ibn == i)[0] for i in xrange(1, xbins.size)]
Z2_idx = np.concatenate([bin[:counts.min()] for bin in idx])
plt.hist(Z2, n, histtype = "stepfilled", color = "#BFBFBF", lw = 0, range = (xi, xf), normed = True)
plt.hist(Z2[Z2_idx], n, histtype = "step", color = "#FFA500", lw = 2, range = (xi, xf), normed = True)
plt.plot([0.2, 2.5], [1/2.3, 1/2.3], "-r", lw = 2)


#mu_pdf = lambda x : st.norm.pdf(x, 0.3, 1.27)*1.27*np.sqrt(2*np.pi)/1.5
mu_pdf = lambda x : st.norm.pdf(x, 0.3, 0.36)*0.36*np.sqrt(2*np.pi)/0.62

idx  = np.concatenate((Z1_idx, Z2_idx))
mu   = mu[idx]
y    = np.random.rand(mu.size) * 1.8
mask = (y <= mu_pdf(mu))

mm = mask & (mu<=1.)
print mm.sum()

plt.figure()
x = np.linspace(0.1, 1., 1000)
plt.hist(mu, 10, normed = True, histtype = "stepfilled", color = "#BFBFBF", lw = 0, range = (0.1, 1))
plt.hist(mu[mm], 10, normed = True, histtype = "step", color = "#FFA500", lw = 2, range = (0.1, 1))
plt.plot(x, mu_pdf(x), "-r", lw = 2)
plt.xlim(0.1, 1.)
plt.xlabel("mu")
plt.show()
