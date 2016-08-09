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

table = np.loadtxt("../../outputs/photometric_fit/remote_set1/phot_fit.physical")
table = binned_stat(table, 100)

ave_sed = np.loadtxt("../../outputs/spectroscopic_fit/ave_table_din.v3.log")

res_pho = [residual(table[:, i], table[:, i + 1], True) if i == 0 or i == 6 else residual(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2) if i != 6]
res_sed = [residual(table[:, ::2][:, i], ave_sed[:, i], True) if i == 0 or i == 3 else residual(table[:, ::2][:, i], ave_sed[:, i], False) for i in xrange(5) if i != 3]

res_pho.pop(2)
res_sed.pop(2)

med_pho = map(np.median, res_pho)
P16_pho = map(st.scoreatpercentile, res_pho, [16]*len(res_pho))
P84_pho = map(st.scoreatpercentile, res_pho, [84]*len(res_pho))

med_sed = map(np.median, res_sed)
P16_sed = map(st.scoreatpercentile, res_sed, [16]*len(res_sed))
P84_sed = map(st.scoreatpercentile, res_sed, [84]*len(res_sed))

med_wui = np.ma.fix_invalid([10**-0.02-1, -0.03, -0.29])
P16_wui = med_wui - np.ma.fix_invalid([10**-0.11-1, -0.14, -0.30])
P84_wui = med_wui + np.ma.fix_invalid([10**0.06-1, 0.12, 0.32])

med_lee = np.ma.fix_invalid([-0.039, 1.075, np.nan])
P16_lee = med_lee - np.ma.fix_invalid([0.232, 0.510, np.nan])
P84_lee = med_lee + np.ma.fix_invalid([0.232, 0.510, np.nan])

#med_pfo = np.ma.fix_invalid([10**-0.51-1, -1.95, 0.34])
#P16_pfo = np.ma.fix_invalid([10**-1.44-1, -2.76, 0.06])
#P84_pfo = np.ma.fix_invalid([10**-0.13-1, -0.28, 0.66])
med_pfo = np.ma.fix_invalid([10**-0.23-1, -0.31, 0.34])
P16_pfo = np.ma.fix_invalid([10**-0.47-1, -0.74, 0.06])
P84_pfo = np.ma.fix_invalid([10**-0.03-1, 0.03, 0.66])

med_mit = np.ma.fix_invalid([10**0.14-1, np.nan, np.nan])
P16_mit = np.ma.fix_invalid([10**0.10-1, np.nan, np.nan])
P84_mit = P16_mit

lab    = [r"$\mathbf{\Delta\,M/M_{\bf SSAG}}$", r"$\mathbf{\Delta\,\left<\log(t)\right>_M}$", r"$\mathbf{\Delta\,A_V}$"]
#labels = ["", r"${\bf full\,SED}$", r"${\bf photo\,SED}$", r"${\bf Wuyts+09}$", r"${\bf Lee+09}$", r"${\bf Pforr+12}$", r"${\bf Mitchell+13}$", ""]
labels = ["", r"${\bf full\,SED}$", r"${\bf photo\,SED}$", r"${\bf Lee+09}$", r"${\bf Pforr+12}$", r"${\bf Mitchell+13}$", ""]

fig, axs = plt.subplots(nrows = len(lab), ncols = 1, figsize = (4, 9), sharex = True)
plt.xlim(0, len(labels) - 1)

x_data = range(1, len(labels) - 1)
for i in xrange(axs.size) :
	axs[i].grid(True, color = "#4D4D4D")
	axs[i].set_ylabel(lab[i], fontsize = 14)
	axs[i].set_ylim(-1.5,+1.5)
	if i > 0 : axs[i].set_yticks(axs[i].get_yticks()[:-1])
	axs[i].axhline(ls = "--", color = "#1A1A1A")

	#y_data = np.array([med_sed[i], med_pho[i], med_wui[i], med_lee[i], med_pfo[i], med_mit[i]])
	#y_erro = np.array([[P16_sed[i], P16_pho[i], P16_wui[i], P16_lee[i], P16_pfo[i], P16_mit[i]],
	                   #[P84_sed[i], P84_pho[i], P84_wui[i], P84_lee[i], P84_pfo[i], P84_mit[i]]])
	y_data = np.array([med_sed[i], med_pho[i], med_lee[i], med_pfo[i], med_mit[i]])
	y_erro = np.array([[P16_sed[i], P16_pho[i], P16_lee[i], P16_pfo[i], P16_mit[i]],
	                   [P84_sed[i], P84_pho[i], P84_lee[i], P84_pfo[i], P84_mit[i]]])
	y_erro = np.abs(y_erro -  y_data)

	axs[i].errorbar(x_data[0], y_data[0], yerr = y_erro[:, 0].reshape((2, 1)), fmt = "x", ecolor = "#7F7F7F", elinewidth = 2, color = "#7F7F7F", mew = 2)
	axs[i].errorbar(x_data[1], y_data[1], yerr = y_erro[:, 1].reshape((2, 1)), fmt = "x", ecolor = "r", elinewidth = 2, color = "r", mew = 2)
	axs[i].errorbar(x_data[2:], y_data[2:], yerr = y_erro[:, 2:], fmt = "x", ecolor = "k", elinewidth = 2, color = "k", mew = 2)

axs[i].set_xticklabels(labels, rotation = 45)
plt.tight_layout()
fig.subplots_adjust(hspace = 0.0)
plt.savefig("../../plots/photometric_fit/we_versus_others.pdf", bbox_inches = "tight")
#plt.savefig("../../plots/photometric_fit/we_versus_others.png", bbox_inches = "tight")
