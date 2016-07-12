#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import scipy.interpolate as inv
import pyfits as pyf
import os

from matplotlib import cm
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def nbins(sample, range_ = None) :
  IQR = lambda x    : st.scoreatpercentile(x, 75.0) - st.scoreatpercentile(x, 25.0)
  if range_ is None : mn, mx = sample.min(), sample.max()
  else              : mn, mx = range_

  mask    = (sample >= mn) & (sample <= mx)
  binsize = (2 * IQR(sample[mask]) / mask.sum() ** (1. / 3))

  return (mx - mn) / binsize, mn, mx, binsize

# READ DATA ========================================================================================

# read SDSS catalog --------------------------------------------------------------------------------

sdss_cat = "../../../TEG/data/photometry/total_photometry.txt"
u_obs, g_obs, r_obs, i_obs, z_obs, redshift = np.loadtxt(sdss_cat, usecols = range(1, 6) + [11], unpack = True)

# read SSAG catalog --------------------------------------------------------------------------------
#logfile = "../../inputs/SFHs_set2/thread0.1.log"
#u_lib, g_lib, r_lib, i_lib, z_lib = np.loadtxt(logfile, usecols = range(14, 19), unpack = True)
#name = np.loadtxt(logfile, usecols = (0,), dtype = "|S")

logfile = "../../inputs/run09_sfh_catalog.txt"
u_lib, g_lib, r_lib, i_lib, z_lib = np.loadtxt(logfile, usecols = range(39, 39+5), unpack = True)
name = np.loadtxt(logfile, usecols = (0,), dtype = "|S")
samn = np.loadtxt("nnew_sample.txt", dtype = "|S")

table = np.loadtxt(logfile, usecols = (5, 12, 16, 17, 18, 25, 38))
#table = np.loadtxt(logfile, usecols = (1, 7, 11, 12, 13, 27))
#table = np.loadtxt(logfile, usecols = (1, 7, 11, 12, 13)) # <---- Plain
# tform_sam, Z_sam, mass_sam, massb_sam, mass1gyr_sam, Av_sam

# define masks -------------------------------------------------------------------------------------

zmask = redshift < 0.03
lmask = (u_lib - r_lib>0.5)&(u_lib - r_lib<3.)&(g_lib - r_lib<0.94)
nmask = (u_lib - r_lib>0.5)&(u_lib - r_lib<3.)&(g_lib - r_lib<0.81)

# MAKE THE PLOTS ===================================================================================

ocolor_x = (u_obs - g_obs)[redshift < 0.03]
ocolor_y = (g_obs - r_obs)[redshift < 0.03]

lcolor_x = (u_lib - g_lib)
lcolor_y = (g_lib - r_lib)

ocolor = (u_obs - r_obs)[zmask]
lcolor = (u_lib - r_lib)[lmask]

#y, x, patches = plt.hist(ocolor, nbins(ocolor, range_ = (0.5, 3.5))[0], range = (0.5, 3.5), color = "k", normed = True)

idx  = np.argsort(lcolor)

blus = st.norm.pdf(lcolor[idx], loc = 1.5, scale = 0.25)
blus = blus / sum(blus)

reds = st.norm.pdf(lcolor[idx], loc = 2.5, scale = 0.3)
reds = reds / sum(reds)

plt.plot(lcolor[idx], reds, ",r", lcolor[idx], blus, ",b")

# define plot properties ---------------------------------------------------------------------------

fig = plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0, 1.25))
#ax1 = plt.subplot(111, xlim = (0.6, 2.1), ylim = (0.1, 1.0))
axins1 = inset_axes(ax1, width = "2.5%", height = "40%", loc = 4)

ax1.set_xlabel(r"$\mathbf{(u-g)\,[mag]}$", fontsize = 15)
ax1.set_ylabel(r"$\mathbf{(g-r)\,[mag]}$", fontsize = 15)

# plot density contours of the SDSS-DR7 galaxies ---------------------------------------------------

nx, xi, xf, bsx = nbins(ocolor_x)
ny, yi, yf, bsy = nbins(ocolor_y)

print("bin sizes:")
print("xbins : {0:5.3f}".format(bsx))
print("ybins : {0:5.3f}".format(bsy))

H, xe, ye = np.histogram2d(ocolor_x, ocolor_y, bins = (nx, ny), normed = True)
H         = gaussian_filter(H, sigma = 2)

cont = ax1.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r, levels = np.arange(0, 6, 0.5))
cb   = plt.colorbar(cont, cax = axins1)
cb.set_label(r"${\bf galaxies\,per\,bin}$", fontsize = 14)

#ax1.scatter(ocolor_x, ocolor_y, s = 1, c = "#7F7F7F", lw = 0)

# plot 0.1 models per bin contour of SSAG ----------------------------------------------------------

nx, xi, xf, bsx = nbins(lcolor_x)
ny, yi, yf, bsy = nbins(lcolor_y)

print("bin sizes:")
print("xbins : {0:5.3f}".format(bsx))
print("ybins : {0:5.3f}".format(bsy))

H, xe, ye = np.histogram2d(lcolor_x, lcolor_y, bins = (nx, ny), normed = True)
H         = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], linewidths = 1.5, colors = ["#1A1A1A"])
#cont = ax1.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r)
#cb   = plt.colorbar(cont, cax = axins1)
#cb.set_label(r"${\bf galaxies\,per\,bin}$", fontsize = 14)

# plot sample model galaxies -----------------------------------------------------------------------

#sb   = np.random.choice(range(lcolor.size), size = 65, p = blus, replace = False)
#sr   = np.random.choice(range(lcolor.size), size = 35, p = reds, replace = False)

#idx  = np.argsort((u_lib - r_lib)[nmask])
#reds = st.norm.pdf((u_lib - r_lib)[nmask][idx], loc = 2.8, scale = 0.1)
#reds = reds / sum(reds)
#sr_new = np.random.choice(range((u_lib - r_lib)[nmask].size), size = 20, p = reds, replace = False)

sb = [i for i in xrange(name[lmask].size) if name[lmask][i] in samn[:65]]
sr = [i for i in xrange(name[lmask].size) if name[lmask][i] in samn[65:]]

sidx = np.concatenate((sb, sr))
#np.savetxt("new_sample.txt", name[lmask][idx][sidx], "%64s")

bcolor_x = lcolor_x[lmask][sb]
bcolor_y = lcolor_y[lmask][sb]

rcolor_x = lcolor_x[lmask][sr]
rcolor_y = lcolor_y[lmask][sr]

#ecolor_x = lcolor_x[nmask][idx][sr_new]
#ecolor_y = lcolor_y[nmask][idx][sr_new]

sizes  = 20 #(mass1gyr_sam / mass_sam) * 100
colors = "none"#(mass1gyr_sam / mass_sam) * 100#Av_sam #tburst_sam * 1e-9

#sc = ax1.scatter(lcolor_x, lcolor_y, c = colors, s = sizes)
#sc = ax1.scatter(bcolor_x, bcolor_y, lw = 0, c = "#154D85", s = sizes)
#sc = ax1.scatter(rcolor_x, rcolor_y, lw = 0, c = "#BC334B", s = sizes)
sc = ax1.scatter(np.append(bcolor_x, rcolor_x), np.append(bcolor_y, rcolor_y), lw = 0, c = "#000080", s = sizes)
#sc = ax1.scatter(ecolor_x, ecolor_y, lw = 0, c = "#800080", s = sizes)
#ax1.scatter((u_lib - g_lib)[bmask], (g_lib - r_lib)[bmask], lw = 0, c = "k")
#cb = plt.colorbar(sc, cax = axins1)
#cb.set_label("tburst[Gyr]")

# plot colorbar ------------------------------------------------------------------------------------

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")
plt.tight_layout()

# Some relevant parameter distributions of the samples ---------------------------------------------

fig, axs = plt.subplots(ncols = 2, nrows = 2)
axs = np.ravel(axs)
#
tform_sam, Z_sam, mass_sam, massb_sam, mass1gyr_sam, pV_sam, V_sam = table[lmask].T
Av_sam = V_sam - pV_sam
#tform_sam, Z_sam, mass_sam, massb_sam, mass1gyr_sam = table[lmask].T
#Av_sam = np.zeros(tform_sam.size)
#tform_sam    = table[:, 0]
#Z_sam        = table[:, 1]
#mass_sam     = table[:, 2]
#massb_sam    = table[:, 3]
#mass1gyr_sam = table[:, 4]
#Av_sam       = table[:, 5]
#
labs = ["visual extinction (mag)", "formation time (Gyr)", "metallicity (solar)", "mass percent formed in the last Gyr"]
pars = [[Av_sam[sb], Av_sam[sr]], [tform_sam[sb]/1e9, tform_sam[sr]/1e9], [Z_sam[sb], Z_sam[sr]], [(mass1gyr_sam/mass_sam*100)[sb], (mass1gyr_sam/mass_sam*100)[sr]]]
#
for i in xrange(4) :
  if i == 3 : axs[i].hist(pars[i], 10, color = ["#154D85", "#BC334B"], histtype = "barstacked", range = (0, 100))
  else      : axs[i].hist(pars[i], 10, color = ["#154D85", "#BC334B"], histtype = "barstacked")
  axs[i].set_xlabel(labs[i], fontsize = 11)
  if i == 0 or i == 2 : axs[i].set_ylabel(r"counts", fontsize = 11)
#
plt.tight_layout()
