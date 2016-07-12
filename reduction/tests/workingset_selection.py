#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as st
import pyfits as pyf
import os

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
u_obs, g_obs, r_obs, i_obs, z_obs, redshift = np.loadtxt(sdss_cat, usecols = range(1, 6) + [11],
                                                         unpack = True)

# read sample properties ---------------------------------------------------------------------------

sfh_name_sam = np.loadtxt("../../inputs/SFHs/fits_list.log", dtype = "|S")
#fits_list = (os.path.join(root, file) for root, subs, files in os.walk("../../run08/.") for file in files if file.endswith(".fits.gz"))

u_sam, g_sam, r_sam, i_sam, z_sam = [], [], [], [], []
tm_sam, tf_sam, tform_sam, Z_sam, Av_sam, mass_sam, mass1gyr_sam, A_sam, tburst_sam = [], [], [], [], [], [], [], [], []
gamma_sam = []
for name in sfh_name_sam :
#for name in fits_list :
  f = pyf.open("../../inputs/SFHs/" + name)
  #f = pyf.open(name)
  
  u_sam.append(f[0].header["umag"])
  g_sam.append(f[0].header["gmag"])
  r_sam.append(f[0].header["rmag"])
  i_sam.append(f[0].header["imag"])
  z_sam.append(f[0].header["zmag"])

  #tm_sam.append(f[0].header["mwla"])
  #tf_sam.append(f[0].header["rfwla"])
  tm_sam.append(f[0].header["mwage"])
  tf_sam.append(f[0].header["rfwage"])
  tform_sam.append(f[0].header["tform"])
  Z_sam.append(f[0].header["z"])
  Av_sam.append(f[0].header["v"] - f[0].header["pv"])
  mass_sam.append(f[0].header["mass"])
  mass1gyr_sam.append(f[0].header["mass1gyr"])
  A_sam.append(f[0].header["a"])
  tburst_sam.append(f[0].header["burstage"])
  gamma_sam.append(f[0].header["gamma"])

u_sam = np.array(u_sam)
g_sam = np.array(g_sam)
r_sam = np.array(r_sam)
i_sam = np.array(i_sam)
z_sam = np.array(z_sam)

tm_sam       = np.array(tm_sam      )
tf_sam       = np.array(tf_sam      )
tform_sam    = np.array(tform_sam   )
Z_sam        = np.array(Z_sam       )
Av_sam       = np.array(Av_sam      )
mass_sam     = np.array(mass_sam    )
A_sam        = np.array(A_sam       )
tburst_sam   = np.array(tburst_sam  )
mass1gyr_sam = np.array(mass1gyr_sam)

# define masks -------------------------------------------------------------------------------------

zmask = redshift < 0.03

SF_mask = (mass1gyr_sam/mass_sam>0.4) & (Av_sam>1.5) & (Z_sam>1)
PS_mask = (mass1gyr_sam/mass_sam<0.4) & (Av_sam<1.5) & (Z_sam<1)

# MAKE PLOT ========================================================================================

# define color-color plane -------------------------------------------------------------------------

ocolor_x = (u_obs - g_obs)[zmask]
ocolor_y = (g_obs - r_obs)[zmask]
#ocolor_y = (r_obs - i_obs)[zmask]

scolor_x = (u_sam - g_sam)
scolor_y = (g_sam - r_sam)
#scolor_y = (r_sam - i_sam)

# define plot properties ---------------------------------------------------------------------------

fig = plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.6, 2.1), ylim = (0.1, 1.0))
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

sizes  = 20
colors = mass1gyr_sam/mass_sam

sc = ax1.scatter(scolor_x, scolor_y, lw = 0, c = colors, s = sizes, alpha = 0.5)
sc = ax1.scatter(scolor_x[SF_mask], scolor_y[SF_mask], lw = 0, c = "#1E90FF", s = sizes)
sc = ax1.scatter(scolor_x[PS_mask], scolor_y[PS_mask], lw = 0, c = "#FF0000", s = sizes)
#cb = plt.colorbar(sc, cax = axins1)
#cb.set_label(r"${\bf mass fraction}$", fontsize = 14)

# plot colorbar ------------------------------------------------------------------------------------

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

# save figure --------------------------------------------------------------------------------------

plt.tight_layout()
#plt.savefig("../../plots/photometric_fit/working_set.pdf", bbox_inches = "tight")
#plt.savefig("../../plots/photometric_fit/working_set.eps", bbox_inches = "tight")
#plt.savefig("../../plots/photometric_fit/working_set.png", bbox_inches = "tight")

# Some relevant parameter distributions of the samples ---------------------------------------------

#fig, axs = plt.subplots(ncols = 2, nrows = 2)
#axs = np.ravel(axs)
#
#labs = ["visual extinction (mag)", "mass", "mass weighted mean age (log)", "mass percent formed in the last Gyr"]
#pars = [[Av_sam, Av[bmask]], [mass_sam, mass[bmask]], [tm_sam, np.log10(tm[bmask])], [mass1gyr_sam/mass_sam*100, (mass1gyr/mass*100)[bmask]]]
#
#from itertools import product
#for i in xrange(4) :
  #axs[i].hist(pars[i], 10, color = ["k", "navy"], histtype = "barstacked")
  #axs[i].set_xlabel(labs[i], fontsize = 11)
  #if i == 0 or i == 2 : axs[i].set_ylabel(r"counts", fontsize = 11)
#
#plt.tight_layout()
#plt.savefig("../../plots/photometric_fit/some_pars_dist.png")
