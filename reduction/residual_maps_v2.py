#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import pyfits as pyf
import thesis_redux_tools as rt
from matplotlib.mlab import griddata
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.ndimage.filters import gaussian_filter
from matplotlib import rc, cm
from matplotlib.colors import ListedColormap as LC
import sys

sdss_cat = "../../../degree_thesis/data/photometry/total_photometry.txt"
u_obs, g_obs, r_obs, i_obs, z_obs, redshift = np.loadtxt(sdss_cat, usecols = range(1, 6) + [11], unpack = True)

name_sed  = np.loadtxt("../../outputs/spectroscopic_fit/names.log", dtype = "|S")
name_pho  = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", usecols = (0,), dtype = "|S")
order_sed = [j for i in xrange(120) for j in xrange(12000) if name_pho[i] == name_sed[j][:12]+".fits.gz"]

table = np.loadtxt("../../outputs/photometric_fit/remote_set3/photofit_SDSS.physical")

ave_sed = rt.binned_stat(np.loadtxt("../../outputs/spectroscopic_fit/table_din.v4.log")[order_sed], bin_size = 1)
chis    = ave_sed[:, 5:]
ave_sed = ave_sed[:, :5]

sfhs_cat = "../../inputs/run09_sfh_catalog.txt"
u_lib, g_lib, r_lib, i_lib, z_lib, tm, tf, tform, Z, V, pV, mass, mass1gyr, A, mu = np.genfromtxt(sfhs_cat,
missing = '""', usecols = range(39, 44) + [33, 35, 5, 12, 38, 25, 16, 18, 19, 14], unpack = True)

u, g, r, i, z = np.genfromtxt("../../inputs/set3_catalog.txt", usecols = range(39, 39 + 5), unpack = True, missing = '""')

table[:, 6:8] = np.log10(table[:, 6:8])

res_sed = [rt.err(table[:, ::2][:, i], ave_sed[:, i]) if i in [0] else rt.err(table[:, ::2][:, i], ave_sed[:, i], False) for i in xrange(5)]

res_sed.pop(2)

lcolor_x = (u_lib - g_lib)[((V - pV) < 1.2) & (mu < 1.)]
lcolor_y = (g_lib - r_lib)[((V - pV) < 1.2) & (mu < 1.)]

zmask = redshift < 0.03
ocolor_x = (u_obs - g_obs)[zmask]
ocolor_y = (g_obs - r_obs)[zmask]

ug = np.repeat(u - g, 100)
gr = np.repeat(g - r, 100)

nx, xi, xf, bsx = rt.nbins(ocolor_x)
ny, yi, yf, bsy = rt.nbins(ocolor_y)

H, xe, ye = np.histogram2d(ocolor_x, ocolor_y, bins = (nx, ny), normed = True)
H         = gaussian_filter(H, sigma = 2)

lab = [r"$M/M_\textrm{SSAG}$", r"$\left<\log(t)\right>_M$", r"$\left<\log(Z/Z_\odot)\right>_M$", r"$A_V$"]


plt.figure(figsize = (9, 8))
axs = [plt.subplot2grid((16, 17), (0, 0), rowspan = 8, colspan = 8),
       plt.subplot2grid((16, 17), (0, 8), rowspan = 8, colspan = 8),
       plt.subplot2grid((16, 17), (8, 0), rowspan = 8, colspan = 8),
       plt.subplot2grid((16, 17), (8, 8), rowspan = 8, colspan = 8)]

axc = plt.subplot2grid((16, 17), (0, 16), rowspan = 16, colspan = 1)
axins1 = inset_axes(axs[0], width = "2.5%", height = "40%", loc = 4)


colx = np.arange(ug.min(), ug.max(), 0.002)
coly = np.arange(gr.min(), gr.max(), 0.002)

nx, xi, xf, bsx = rt.nbins(lcolor_x)
ny, yi, yf, bsy = rt.nbins(lcolor_y)

#6600FF
#0066CC

H2, xe2, ye2 = np.histogram2d(lcolor_x, lcolor_y, bins = (nx, ny), normed = True)
H2           = gaussian_filter(H2, sigma = 2)

blues = LC(np.loadtxt("my_blues.dat")/255.)
krp   = LC(np.loadtxt("kriptonite.dat")/255.)
for i in xrange(len(axs)) :
  axs[i].set_xlim(0.67, 2.05)
  axs[i].set_ylim(0.125, 0.96)
  axs[i].set_title(lab[i])

  cont = axs[i].contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r, levels = np.arange(0, 6, 0.5))

  #r = axs[i].hexbin(ug[::100], gr[::100], rt.binned_stat(res_sed[i], 100), gridsize = 40, cmap = cm.afmhot_r, vmin = -0.3, vmax = 0.3)
  r = axs[i].hexbin(ug, gr, res_sed[i], gridsize = 40, cmap = cm.afmhot_r, vmin = -0.3, vmax = 0.3)

  if i in [0, 1] : axs[i].tick_params(labelbottom = False)
  if i in [1, 3] : axs[i].tick_params(labelleft = False)

  if i in [2, 3] : axs[i].set_xlabel("$u-g$", fontsize = 16)
  if i in [0, 2] : axs[i].set_ylabel("$g-r$", fontsize = 16)

cb   = plt.colorbar(cont, cax = axins1)
cb.set_label("galaxies per bin", fontsize = 14)
axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

cb = plt.colorbar(r, cax = axc)
cb.set_label(r"Median $\Delta$", fontsize = 16)

plt.tight_layout()
plt.savefig("residual_maps.eps", bbox_inches = "tight")
plt.show()
