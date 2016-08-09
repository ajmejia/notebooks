#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import DataLoader as dl
import scipy.stats as st
from itertools import product
from matplotlib import rc

rc("lines", linewidth=1.0)
rc("patch", linewidth=1.0)
rc("font", family="serif", serif="Times New Roman", size=9.0)
rc("savefig", dpi=92)
rc("axes", linewidth=0.5, labelsize=9.0, titlesize=9.0)
rc("legend", fontsize="small")
rc("xtick.major", width=0.3)
rc("xtick", labelsize="small")
rc("ytick.major", width=0.3)
rc("ytick", labelsize="small")

z      = np.loadtxt("../../inputs/universe_age.dat", dtype = np.str, usecols = (0,))
age    = np.loadtxt("../../inputs/universe_age.dat", usecols = (1,))
IDs    = np.loadtxt("../../inputs/photoz3/set3.counts", dtype = np.str, usecols = (0,))
counts = np.loadtxt("../../inputs/photoz3/set3.counts", dtype = np.int, usecols = (1,))
output = ["z0p50", "z1p00", "z1p50", "z2p00", "z2p50"]
#output = ["z1p50", "z2p00", "z2p50", "z3p00"]
mask   = np.array([True if ID.replace(".", "p") in output else False for ID in IDs], dtype = np.bool)
data_z = [dl.load_data(ID.replace(".", "p"), count, 50, 56) for ID, count in zip(IDs[mask], counts[mask])]
ages   = [age[i] for i in xrange(z.size) if "z" + z[i].replace(".", "p") in output]
lab    = [r"$\Delta M_\star/M_\star^\text{SSAG}$", r"$\Delta\,\left<\log{t_\star}\right>_M$", r"$\Delta\,\left<\log{t_\star}\right>_L$", r"$\Delta\,\left<\log{Z_\star/Z\odot}\right>_M$", r"$\Delta\,A_V$"]
lb     = "M log_t_M log_t_L log_Z_M Av".split()

fig, axs = plt.subplots(len(output), 5, sharex = True, figsize = (7, 5))

for i, j in product(xrange(len(output)), xrange(5)) :

  mask = data_z[i].physical["log_t_M_mod"] < np.log10(ages[i] * 1e9)

  med = np.median(data_z[i].residuals[lb[j]][mask])
  p16 = st.scoreatpercentile(data_z[i].residuals[lb[j]][mask], 16)
  p84 = st.scoreatpercentile(data_z[i].residuals[lb[j]][mask], 84)

  axs[i, j].hist(data_z[i].residuals[lb[j]][mask], 40, histtype = "stepfilled", alpha = 0.5, ec = "#0062FF", fc = "#0062FF", range = (-1.5, +1.5), lw = 2)
  axs[i, j].axvline(med, ls = "--", color = "#2e3436")

  st_textp = r"$%+.2f_{%+.2f}^{%+.2f}$" % (med, p16, p84)
  axs[i, j].text(0.02, 0.70, st_textp, fontsize="x-small", color="k", ha="left", transform=axs[i, j].transAxes)

  if i == len(output) - 1 : axs[i, j].set_xlabel(lab[j])
  if j == 0               : axs[i, j].set_ylabel(r"${0}={1}$".format(output[i][0], output[i][1:].replace("p", ".")))

  axs[i, j].set_xlim(-1.5, +1.5)
  axs[i, j].set_xticklabels(["", -1.0, "", 0.0, "", +1.0, ""])
  axs[i, j].set_yticks([])

plt.tight_layout()
plt.subplots_adjust(hspace = 0.0)
#plt.savefig("residuals_z.pdf", bbox_inches="tight", facecolor="#C4C4C0", edgecolor="none")
plt.savefig("residuals_z.eps", bbox_inches="tight")
#plt.show()
