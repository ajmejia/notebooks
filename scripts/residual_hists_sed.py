#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import pyfits as pyf
import thesis_redux_tools as rt

def binned_stat(table, bin_size, stat = "mean") :
	'''
	This function computes the given statistic in a table (row sense)
	by bins of a given size.
	
	input parameters:
	----------------
	*  table    : an arraylike table to be binned_averaged.
	*  bin_size : the number of elements in each bin. It's mandatory
	              that bin_size * int(# of bins) = table.shape[0].
  *  stat     : character string with the name of the statistic to compute. Available options are:
                - mean (default)
                - median
                - stdev
	'''
	N_bins = table.shape[0] / bin_size
	if N_bins * bin_size == table.shape[0] :
		if stat == "mean" :
			return np.array([np.mean(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
		elif stat == "median" :
			return np.array([np.median(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
		elif stat == "stdev" :
			return np.array([np.std(table[i*bin_size:(i+1)*bin_size], 0) for i in xrange(N_bins)])
	else :
		raise ValueError, "bin_size must fulfill: bin_size * int(# of bins) = table.shape[0]."

def residual(true_val, retr_val, relative = True) :
	if relative : return (retr_val - true_val) / true_val
	else        : return retr_val - true_val

def Q(true_val, retr_val) :
	return np.sqrt(np.sum((retr_val - true_val) ** 2))

table = np.loadtxt("../../outputs/photometric_fit/remote_set1/phot_fit.physical")
table = binned_stat(table, 100)

table[:, 6:8] = np.log10(table[:, 6:8])

#mass_sam, mass_mod, tm_sam, tm_mod, tf_sam, tf_mod, Z_sam, Z_mod, Av_sam, Av_mod = table.T
ave_sed = np.loadtxt("../../outputs/spectroscopic_fit/ave_table_din.v3.log")

res_pho = [residual(table[:, i], table[:, i + 1]) if i == 0 else residual(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2) if i != 4]
res_sed = [residual(table[:, ::2][:, i], ave_sed[:, i]) if i == 0 else residual(table[:, ::2][:, i], ave_sed[:, i], False) for i in xrange(5) if i != 2]

avr_pho = map(np.mean, res_pho)
P16_pho = map(st.scoreatpercentile, res_pho, [16]*len(res_pho))
P84_pho = map(st.scoreatpercentile, res_pho, [84]*len(res_pho))
std_pho = map(np.std, res_pho)

avr_sed = map(np.mean, res_sed)
P16_sed = map(st.scoreatpercentile, res_sed, [16]*len(res_sed))
P84_sed = map(st.scoreatpercentile, res_sed, [84]*len(res_sed))
std_sed = map(np.std, res_sed)

#lab = ["mass residuals", "mass weighted mean age residuals", "flux weighted mean age residuals", "metallicity residuals", "visual extinction residuals"]
#lab = ["mass residuals", "mass weighted mean age residuals", "flux weighted mean age residuals", "visual extinction residuals"]
lab = [r"$\mathbf{\Delta\,M/M_{\bf SSAG}}$", r"$\mathbf{\Delta\,\left<\log(t)\right>_M}$", r"$\mathbf{\Delta\,\left<\log(Z/Z_\odot)\right>_M}}$", r"$\mathbf{\Delta\,A_V}$"]

print "{0:^47s}".format("photometric fit")
for i in xrange(len(avr_pho)) : print "{0:32s} {1:5.2f} {2:5.2f} {3:5.2f}".format(lab[i], avr_pho[i], P16_pho[i], P84_pho[i])
print
print "{0:^47s}".format("spectroscopic fit")
for i in xrange(len(avr_sed)) : print "{0:32s} {1:5.2f} {2:5.2f} {3:5.2f}".format(lab[i], avr_sed[i], P16_sed[i], P84_sed[i])

fig, axs = plt.subplots(ncols = 2, nrows = 2, sharex = "col", figsize = (8, 7))
axs = np.ravel(axs)

for i in xrange(4) :
	ran = axs[i].set_xlim(-1.5, +1.5)
	axs[i].set_ylim(0,50)
	axs[i].set_xlabel(lab[i], fontsize = 14)

	axs[i].text(-1.4, 45, "$\mu={0:5.2f}$".format(avr_pho[i]), fontsize = 14, color = "r")
	axs[i].text(-1.4, 39, "$\mu={0:5.2f}$".format(avr_sed[i]), fontsize = 14)
	#axs[i].text(-1.4, 46, "$\mu={0:5.2f}$".format(avr_sed[i]), fontsize = 14)
	#axs[i].text(-1.4, 40, "$\sigma={0:5.2f}$".format(std_sed[i]), fontsize = 14)
	fig.text(0.05, 0.55, r"$\bf Frequency$", fontsize = 14, rotation = "vertical")
	
	#axs[i].set_yticks([])
	if i == 2 : axs[i].set_xticks(axs[i].get_xticks()[:- 1])
	if i == 1 or i == 3 : axs[i].tick_params(labelleft = False, labelright = True)

	if i == 2 :
		axs[i].hist(res_sed[i], 10, ec = "#7F7F7F", fc = "#7F7F7F")
		axs[i].hist(res_pho[i], 10, ec = "k", fc = "none", ls = "dotted", histtype = "step")
	else :
		axs[i].hist(res_sed[i], rt.nbins(res_pho[i], range_ = ran)[0], ec = "#7F7F7F", fc = "#7F7F7F")
		axs[i].hist(res_pho[i], rt.nbins(res_pho[i], range_ = ran)[0], ec = "k", fc = "none", ls = "dotted", histtype = "step")
	axs[i].axvline(avr_pho[i], ls = "--", color = "r", lw = 2, label = r"${\rm photo\,SED}$")
	axs[i].axvline(avr_sed[i], ls = "--", color = "k", lw = 2, label = r"${\rm full\,SED}$")
	#axs[i].text(-1.49, 1.0, r"$68\%\,{\rm range}\,{\rm ratio}$", fontsize = 12)

	#if i == 0 : axs[i].legend(loc = "upper left", frameon = False, fontsize = 12)

#plt.tight_layout()
fig.subplots_adjust(wspace = 0.0)
plt.savefig("../../plots/photometric_fit/residual_hists_sed.pdf", bbox_inches = "tight")
#plt.savefig("../../plots/photometric_fit/residual_hists.png", bbox_inches = "tight")
#plt.show()
