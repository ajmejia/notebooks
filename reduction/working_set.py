#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as st
import pyfits as pyf

from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from thesis_redux_tools import nbins

# READ DATA ========================================================================================

# read color redshift evolution --------------------------------------------------------------------

rs, ug_zev_old = np.loadtxt("../../inputs/ug_zev_13gyr", usecols = (0,8,), unpack = True)
gr_zev_old = np.loadtxt("../../inputs/gr_zev_13gyr", usecols = (8,))
ug_zev_new = np.loadtxt("../../inputs/ug_zev_001gyr", usecols = (8,))
gr_zev_new = np.loadtxt("../../inputs/gr_zev_001gyr", usecols = (8,))

# read SDSS catalog --------------------------------------------------------------------------------
sdss_cat = "../../../TEG/data/photometry/total_photometry.txt"
u_obs, g_obs, r_obs, i_obs, z_obs, redshift = np.loadtxt(sdss_cat, usecols = range(1, 6) + [11], unpack = True)

# read sample properties ---------------------------------------------------------------------------

#sfh_name_sam = np.loadtxt("../../inputs/SFHs_set1/set1.log", dtype = "|S")
sfh_name_sam = np.loadtxt("../../inputs/SFHs_set2/fits_list.log", dtype = "|S")

u_sam, g_sam, r_sam, i_sam, z_sam = [], [], [], [], []
tm_sam, tf_sam, tform_sam, Z_sam, Av_sam, mass_sam, mass1gyr_sam, A_sam, tburst_sam = [], [], [], [], [], [], [], [], []
for name in sfh_name_sam :
  f = pyf.open("../../inputs/SFHs_set2/" + name)
  
  u_sam.append(f[0].header["umag"])
  g_sam.append(f[0].header["gmag"])
  r_sam.append(f[0].header["rmag"])
  i_sam.append(f[0].header["imag"])
  z_sam.append(f[0].header["zmag"])

  #tm_sam.append(f[0].header["mwla"])
  #tf_sam.append(f[0].header["rfwla"])
  tform_sam.append(f[0].header["tform"])
  Z_sam.append(f[0].header["z"])
  Av_sam.append(f[0].header["v"] - f[0].header["pv"])
  mass_sam.append(f[0].header["mass"])
  mass1gyr_sam.append(f[0].header["mass1gyr"])
  A_sam.append(f[0].header["a"])
  tburst_sam.append(f[0].header["burstage"])

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

# read SFH library data ----------------------------------------------------------------------------

#blues = np.array(['56421.152274', '56426.486358', '56426.230352', '56427.500769',
       #'56422.690793', '56423.228176', '56426.929471', '56421.140615',
       #'56427.732118', '56428.181762', '56421.550003', '56426.589616',
       #'56428.640047', '56427.990432', '56423.493704', '56424.823774',
       #'56425.924211', '56430.897102', '56430.039940', '56422.653237',
       #'56423.992303', '56423.867777', '56425.149235', '56421.874639',
       #'56422.508166', '56427.544719', '56431.107434', '56422.436490',
       #'56420.275014', '56420.129787'], 
      #dtype='|S12')

blues = np.loadtxt("../../inputs/blues.log", dtype = "|S")

sfhs_cat = "../../inputs/run07_sfh_catalog.txt"
sfh_name_lib = np.genfromtxt(sfhs_cat, missing = '""', usecols = [0], dtype = "|S")
tburst, bext = np.genfromtxt(sfhs_cat, missing = '""', usecols = [6, 7], unpack = True)
u_lib, g_lib, r_lib, i_lib, z_lib, tm, tf, tform, Z, V, pV, mass, mass1gyr, A, mu = np.genfromtxt(sfhs_cat,
missing = '""', usecols = range(39, 44) + [33, 35, 5, 12, 38, 25, 16, 18, 19, 14], unpack = True)
#missing = '""', usecols = range(27, 32) + [33, 35, 5, 12, 38, 25, 16, 18, 19, 14], unpack = True)  # plain fluxes

Av = V - pV

# define masks -------------------------------------------------------------------------------------

zmask = redshift < 0.03
dmask = (Av < 2.0) & (mu <= 1.0)
bmask = np.zeros(r_lib.size, dtype = np.bool)
bmask[np.array([np.where(sfh_name_lib == b)[0][0] for b in blues])] = True

# MAKE PLOT ========================================================================================

# define color-color plane -------------------------------------------------------------------------

ocolor_x = (u_obs - g_obs)[zmask]
ocolor_y = (g_obs - r_obs)[zmask]
#ocolor_y = (r_obs - i_obs)[zmask]

lcolor_x = (u_lib - g_lib)[dmask]
lcolor_y = (g_lib - r_lib)[dmask]
#lcolor_y = (r_lib - i_lib)[dmask]

scolor_x = (u_sam - g_sam)
scolor_y = (g_sam - r_sam)
#scolor_y = (r_sam - i_sam)

# define plot properties ---------------------------------------------------------------------------

fig = plt.figure(figsize = (7, 6))
#ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0, 1.25))
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

#ax1.scatter(ocolor_x, ocolor_y, s = 1, c = "#7F7F7F", lw = 0)

# plot 0.1 models per bin contour of SSAG ----------------------------------------------------------

#nx, xi, xf, bsx = nbins(lcolor_x)
#ny, yi, yf, bsy = nbins(lcolor_y)
#
#print("bin sizes:")
#print("xbins : {0:5.3f}".format(bsx))
#print("ybins : {0:5.3f}".format(bsy))
#
#H, xe, ye = np.histogram2d(lcolor_x, lcolor_y, bins = (nx, ny), normed = True)
#H         = gaussian_filter(H, sigma = 2)
#ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.1], linewidths = 1.5, colors = ["#1A1A1A"])
#cont = ax1.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r)
#cb   = plt.colorbar(cont, cax = axins1)
#cb.set_label(r"${\bf galaxies\,per\,bin}$", fontsize = 14)

# plot sample model galaxies -----------------------------------------------------------------------

#sss    = np.ones(rs.size)*20
#sss[5] = 70

#ax1.scatter(ug_zev_old, gr_zev_old, c = rs, s = sss, lw = 0)
#ax1.scatter(ug_zev_new, gr_zev_new, c = rs, s = sss, lw = 0, marker = "s")

sizes  = 20 #(mass1gyr_sam / mass_sam) * 100
colors = "k"#(mass1gyr_sam / mass_sam) * 100#Av_sam #tburst_sam * 1e-9

sc = ax1.scatter(scolor_x, scolor_y, lw = 0, c = colors, s = sizes)
#ax1.scatter((u_lib - g_lib)[bmask], (g_lib - r_lib)[bmask], lw = 0, c = "k")
#cb = plt.colorbar(sc, cax = axins1)
#cb.set_label("tburst[Gyr]")

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
