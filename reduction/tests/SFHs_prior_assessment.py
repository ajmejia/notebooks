#!/usr/bin/python
#
# Evaluation of the prior distributions for the SAM of galaxy formation, as presented in
# Chen et al. (2012).
#
# tform    : onset of star formation. Uniformly distributed over the range 13-1.5 Gyr.
# gamma    : inverse time scale. Uniformly distributed over the range 0-1 1/Gyr.
# A        : amplitude of the burst (=M_burst/M_continuum). Logaritmically distributed over the
#            range 0.03-4.
# burstext : duration of the burst. Uniformly distributed over the range 30-300 Myr.
# tcut     : time of truncation of the SFR.
# taucut   : time scale of the SFR(t > tcut). The log10(taucut) is uniformly distributed over the
#            range 7-9.
# Z        : metallicity. 95% of the galaxies are uniformly distributed over 0.2-2.5 Zo and the 5%
#            are uniformly distributed over 0.02-0.2 Zo.
# tau      : V-band optical depth. Has a Gaussian distribution over 0-6 with peak at 1.2 and 68% of
#            probability distribution sampled over the range 0-2.
# mu       : fraction of the optical depth affecting populations older than 10 Myr. Has a Gaussian
#            distribution with peak at 0.3 and 68% of the probability distributed over 0.1-1.
# sigma    : line of sight velocity dispersion. Uniformly distributed over 50-400 Km/s.
#
# The SFR is such that:
#   * 15% of the galaxies experiences a burst in the last 2 Gyr.
#   * 30% of the galaxies experiences a truncation at a random time in the past (tcut).
#   * 10% of the galaxies experiences a burst in the last 2 Gyr (after truncations are included).
#
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st

IQR = lambda x        : st.scoreatpercentile(x, 75.0) - st.scoreatpercentile(x, 25.0)

def nbins(sample, range_ = None) :
	if range_ is None : mn, mx = sample.min(), sample.max()
	else              : mn, mx = range_

	mask    = (sample >= mn) & (sample <= mx)
	binsize = (2 * IQR(sample[mask]) / mask.sum() ** (1. / 3))

	return (mx - mn) / binsize, mn, mx, binsize

#tform, gamma, A, burstext, tcut, taucut, Z, tau, mu, sigma = np.genfromtxt("../inputs/run07_sfh_catalog.txt", missing = '""', usecols = (5, 8, 19, 7, 10, 11, 12, 13, 14, 15), unpack = True)
table = np.genfromtxt("../inputs/run09_sfh_catalog.txt", missing = '""', usecols = (5, 8, 19, 7, 10, 11, 12, 13, 14, 15))
trun = np.genfromtxt("../inputs/run09_sfh_catalog.txt", missing = '""', usecols = [9], dtype = np.bool)
tburst = np.genfromtxt("../inputs/run09_sfh_catalog.txt", missing = '""', usecols = [6])

#Z = table[:, 6]
#
#Z1 = Z[Z<=0.2]
#Z2 = Z[Z>0.2]
#table_Z1 = table[Z <= 0.2]
#table_Z2 = table[Z > 0.2]
#trun_Z1 = trun[Z <= 0.2]
#trun_Z2 = trun[Z > 0.2]
#tburst_Z1 = trun[Z <= 0.2]
#tburst_Z2 = trun[Z > 0.2]
#
#n, xi, xf, bs = nbins(Z1, range_ = (0.058920444050700002, 0.2))
#counts, xbins = np.histogram(Z1, n, range = (xi, xf))
#ibn = np.digitize(Z1, xbins)
#idx = [np.where(ibn == i)[0] for i in xrange(1, xbins.size)]
#Z1_idx = np.concatenate([bin[:counts.min()] for bin in idx])
#
#n, xi, xf, bs = nbins(Z2, range_ = (0.2, 2.5))
#counts, xbins = np.histogram(Z2, n, range = (xi, xf))
#ibn = np.digitize(Z2, xbins)
#idx = [np.where(ibn == i)[0] for i in xrange(1, xbins.size)]
#Z2_idx = np.concatenate([bin[:counts.min()] for bin in idx])
#
#trun = np.concatenate((trun_Z1, trun_Z2))
#tburst = np.concatenate((tburst_Z1, tburst_Z2))
#new_table = np.concatenate((table_Z1[Z1_idx, :], table_Z2[Z2_idx, :]))

tform    = table[:, 0]
gamma    = table[:, 1]
A        = table[:, 2]
burstext = table[:, 3]
tcut     = table[:, 4]
taucut   = table[:, 5]
Z        = table[:, 6]
tau      = table[:, 7]
mu       = table[:, 8]
sigma    = table[:, 9]

print table.shape, new_table.shape

p1 = int(round(np.count_nonzero(tburst<2e9)/1e5*100))
p2 = int(round(trun.sum()/1e5*100))
p3 = int(round(np.count_nonzero(Z < 0.2)/1e5*100))
p4 = int(round(np.count_nonzero(Z > 0.2)/1e5*100))

print("galaxies with burst within the last 2 Gyr: ~{0} %".format(p1))
print("galaxies with SFR truncated              : ~{0} %".format(p2))
print
print("galaxies with Z < 0.2                    : {0} %".format(p3))
print("galaxies with Z > 0.2                    : {0} %".format(p4))

#mu_pdf = lambda x : st.norm.pdf(x, 0.3, 1.27)*1.27*np.sqrt(2*np.pi)/1.5
mu_pdf = lambda x : st.norm.pdf(x, 0.3, 0.36)*0.36*np.sqrt(2*np.pi)/0.62

y    = np.random.rand(mu.size) * 1.8
mask = (y <= mu_pdf(mu)) & (mu<=1.0)

tform = tform[mask]
gamma = gamma[mask]
A     = A[mask]
burstext = burstext[mask]
tcut = tcut[mask]
taucut = taucut[mask]
Z = Z[mask]
tau = tau[mask]
mu = mu[mask]
sigma = sigma[mask]

tcut   = tcut[~np.ma.fix_invalid(tcut).mask]
taucut = taucut[~np.ma.fix_invalid(taucut).mask]

pars = [tform/1e9, gamma*1e9, np.log10(A), burstext/1e8, tcut/1e9, np.log10(taucut), Z, tau, mu, sigma]
labs = ["tform (Gyr)", "gamma (1/Gyr)", "log(A)", "burstext (100 Myr)", "tcut (Gyr)", "log(taucut)", "Z (solar)", "tau", "mu", "sigma"]

xs = [np.linspace(1.5, 13.5, 1000), np.linspace(0, 1, 1000), np.linspace(-1.52, 0.60, 1000),
      np.linspace(0.3, 3.0, 1000), np.linspace(1.5, 13.5, 1000), np.linspace(7, 9, 1000),
      np.linspace(0.02, 2.5, 1000), np.linspace(0, 6, 1000), np.linspace(0.1, 1, 1000),
      np.linspace(50, 400, 1000)]
ys = [lambda x : [1./12]*1000, lambda x : [1.]*1000, lambda x : [1./2.12]*1000, lambda x : [1./2.7,]*1000,
      lambda x : [1./12]*1000, lambda x : [1./2]*1000, lambda x : np.piecewise(x, [x < 0.2, x > 0.2], [1./0.14, 1./2.3]),
      lambda x : st.norm.pdf(x, 1.2, 0.98520478)*0.98520478*np.sqrt(2*np.pi)/2.19, lambda x : st.norm.pdf(x, 0.3, 0.36566883)*0.36566883*np.sqrt(2*np.pi)/0.62, lambda x : [1./350]*1000]

fig = plt.figure(figsize = (6, 9))

axs = [fig.add_subplot(5, 2, i + 1) for i in xrange(10)]

for i, ax in enumerate(axs) :
	ax.tick_params(labelleft = False, labelsize = 9)
	range_ = ax.set_xlim(xs[i].min(), xs[i].max())
	ax.set_xlabel(labs[i], fontsize = 9)
	if i == 4 : ax.set_ylabel("pdf", fontsize = 9)

	if i == 6 :
		ax.hist(Z[Z<=0.2], 1, fc = "#FFA500", ec = "none", normed = True)
		ax.hist(Z[Z>0.2], 16, fc = "#FFA500", ec = "none", normed = True)
		ax.hist(Z[Z<=0.2], 1, fc = "none", ec = "#1A1A1A", lw = 2, normed = True)
		ax.hist(Z[Z>0.2], 16, fc = "none", ec = "#1A1A1A", lw = 2, normed = True)
	else :
		if i == 4 or i == 5 :
			ax.hist(pars[i], 10, fc = "#FFA500", ec = "none", normed = True, range = range_)
		else :
			ax.hist(pars[i], 10, fc = "#FFA500", ec = "none", normed = True, range = range_)
		ax.hist(pars[i], 10, fc = "none", ec = "#1A1A1A", lw = 2, normed = True, range = range_)
	func = ys[i]
	ax.plot(xs[i], func(xs[i]), "-r", lw = 2)

fig.tight_layout()
plt.show()
