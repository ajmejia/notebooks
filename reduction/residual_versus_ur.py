#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats.mstats as st
import pyfits as pyf
import thesis_redux_tools as rt
import sys

fid = "SDSS"

#if len(sys.argv) > 1 : fid = sys.argv[1]
#else :
  #print "Usage: residual_versus_ur.py physical_file_id"
  #sys.exit(1)

fits_list = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", dtype = np.str, usecols = (0,))
table     = np.loadtxt("../../outputs/photometric_fit/remote_set3/photofit_" + fid + ".physical")

#table = rt.binned_stat(table, 100)

table = table[table[:, 6] > 0.4, :]

u, g, r, i, z = [], [], [], [], []
sfr = []
for name in fits_list :
  f = pyf.open("../../inputs/SFHs_set3/" + name)

  u.append(f[0].header["umag"])
  g.append(f[0].header["gmag"])
  r.append(f[0].header["rmag"])
  i.append(f[0].header["imag"])
  z.append(f[0].header["zmag"])

  sfr.append(f[0].header["mass1gyr"] / f[0].header["mass"])

u  = np.repeat(u, 100)
g  = np.repeat(g, 100)
r  = np.repeat(r, 100)
i  = np.repeat(i, 100)
z  = np.repeat(z, 100)

sfr = np.array(sfr)

ur = (u - r)[table[:, 6] > 0.4]

table = rt.binned_stat(table, 100, "median")
ur    = rt.binned_stat(ur, 100, "median")

table[:, 6:8] = np.log10(table[:, 6:8])
res = [rt.err(table[:, i], table[:, i + 1]) if i == 0 else rt.err(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2)]
lab = [r"$\Delta\,M/M_\textrm{SSAG}$",
       r"$\Delta\,\left<\log(t)\right>_M$",
       r"$\Delta\,\left<\log(t)\right>_{L_r}$",
       r"$\Delta\,\left<\log(Z/Z_\odot)\right>_M$",
       r"$\Delta\,A_V$"]

fig, axs = plt.subplots(1, len(res), figsize = (20, 3))

mask = (res[0] > 0)
print min(table[:, 2][mask]), max(table[:, 2][mask])
print table[res[1][mask] == res[1][mask].min(), ::2], res[1][mask].min()

for i in xrange(len(res)) :
  axs[i].set_title(lab[i])
  axs[i].set_xlim(0.7, 3.1)
  axs[i].set_ylim(-1.5, +1.5)
  #axs[i].set_ylim(-0.7, +0.7)

  #x, y, yerr = rt.binner(ur, res[i], "median", 10, ebar = True)

  sc = axs[i].scatter(ur[mask], res[i][mask], c = table[mask, 6], lw = 0, s = 15, vmin = table[:, 6].min(), vmax = table[:, 6].max())
  #axs[i].plot(ur, res[i], ".", color = "#7F7F7F", alpha = 0.2)
  #axs[i].plot(x, y - yerr[0], "--k", lw = 1.5)
  #axs[i].plot(x, y + yerr[1], "--k", lw = 1.5)
  #axs[i].plot(x, y, "-k", lw = 1.5)
  axs[i].axhline(ls = "--", color = "#A52A2A", lw = 1.7)

#plt.colorbar(sc, cax = axs[-1])
axs[0].set_ylabel("residuals", fontsize = 14)
axs[2].set_xlabel("$(u-r)$[mag]", fontsize = 14)

#plt.tight_layout()
fig.subplots_adjust(wspace = 0.25, hspace = 0., left = 0.06, right = 0.96, top = 0.85)
#plt.savefig("../../plots/photometric_fit/residual_versus_ur_" + fid + ".png", bbox_inches = "tight")
plt.show()
