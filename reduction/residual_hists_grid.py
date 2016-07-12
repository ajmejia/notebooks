#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats.mstats as st
import pyfits as pyf
import thesis_redux_tools as rt

from operator import __getslice__

#table = np.loadtxt("../../outputs/photometric_fit/remote_set3/phot_fit_SDSS.physical")
#table = np.loadtxt("../../outputs/photometric_fit/phot_fit_SDSS_dust_fitting_monoZ_dust_fitting_off.physical")
table = np.loadtxt("../../outputs/photometric_fit/remote_set3/photofit_SDSS.physical")
names = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", dtype = np.str)
names = map(__getslice__, names, [0] * names.size, [12] * names.size)

table[:, 6:8] = np.log10(table[:, 6:8])

res = [rt.err(table[:, i], table[:, i + 1]) if i == 0 else rt.err(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2)]
lab = ["mass residuals", "mass weighted age residuals", "flux weighted age residuals", "metallicity residuals", "dust extinction residuals"]

for j in xrange(len(res)) :
  fig, axs = plt.subplots(5, 13, sharex = True, sharey = True, figsize = (20, 15))

  plt.xlim(-1.5, +1.5)
  plt.ylim(0, 40)

  axs = np.ravel(axs)

  for i in xrange(65) :
    data   = res[j][i * 100:(i + 1) * 100]
    median = np.median(data)
    p16    = st.scoreatpercentile(data, 16.0)
    p84    = st.scoreatpercentile(data, 84.0)

    counts, bins, patches = axs[i].hist(data, 30, histtype = "step", hatch = "///", lw = 1, color = "#1A1A1A", range = (-1.5, +1.5))

    axs[i].axvline(median, ls = "--", lw = 1.5, color = "#000080")
    axs[i].axvline(p16, ls = "-.", lw = 1.5, color = "#000080")
    axs[i].axvline(p84, ls = "-.", lw = 1.5, color = "#000080")

  axs[-1].set_xticks([-1., 0, +1.])
  axs[52].set_yticks(list(axs[i].get_yticks()[1:-1]))
  axs[58].set_xlabel(lab[j], fontsize = 16)
  axs[26].set_ylabel("counts", fontsize = 16)

  plt.tight_layout()
  plt.subplots_adjust(wspace = 0.01, hspace = 0.01, bottom = 0.06)
  #plt.savefig("../../plots/photometric_fit/" + lab[j].replace(" ", "_") + ".png", bbox_inches = "tight")

plt.show()
