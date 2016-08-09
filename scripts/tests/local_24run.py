#!/usr/bin/python

from matplotlib import cm
import matplotlib.pyplot as plt
import numpy as np
import pyfits as pyf
from scipy.stats.mstats import scoreatpercentile as percentile
from scipy.ndimage.filters import gaussian_filter
from operator import add
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

err = lambda obs, ret : (ret - obs) / obs * 100.0

IQR = lambda x        : percentile(x, 75.0) - percentile(x, 25.0)

def nbins(sample, range_ = None) :
	if range_ is None : mn, mx = sample.min(), sample.max()
	else              : mn, mx = range_

	mask    = (sample >= mn) & (sample <= mx)
	binsize = (2 * IQR(sample[mask]) / mask.sum() ** (1. / 3))

	return (mx - mn) / binsize, mn, mx, binsize

#with open("../inputs/SFHs/fits_list.log", "r") as fits_names : names = [fits_names.readline()[:- 1] for i in xrange(7692)]
#u, g, r, i, z, av, met, mass = [], [], [], [], [], [], [], []
#for name in names :
	#fits_file = pyf.open("../inputs/SFHs/" + name)
#
	#u.append(fits_file[0].header["umag"])
	#g.append(fits_file[0].header["gmag"])
	#r.append(fits_file[0].header["rmag"])
	#i.append(fits_file[0].header["imag"])
	#z.append(fits_file[0].header["zmag"])
#
	#av.append(fits_file[0].header["v"] - fits_file[0].header["pv"])
	#met.append(fits_file[0].header["z"])
	#mass.append(fits_file[0].header["mass"])
#
#u    = np.array(u)
#g    = np.array(g)
#r    = np.array(r)
#i    = np.array(i)
#z    = np.array(z)
#av   = np.array(av)
#met  = np.array(met)
#mass = np.array(mass)
table = np.genfromtxt("../inputs/run07_sfh_catalog.txt", missing = '""', usecols = range(39, 42) + [12, 16, 0, 25, 38])
av    = table[:, - 1] - table[:, - 2]
mask  = av < 2.0
av    = av[mask]
u     = table[mask, 0]
g     = table[mask, 1]
r     = table[mask, 2]
i     = table[mask, 3]
z     = table[mask, 4]
met   = table[mask, 5]
mass  = table[mask, 6]

fitsnames = np.loadtxt("../inputs/SFHs_set1/set1.log", dtype = "|S")

ug_color, gr_color, sav, smet, smass = [], [], [], [], []
for name in fitsnames :
	
	f = pyf.open("../inputs/SFHs_set1/" + name)
	
	ug_color.append(f[0].header["umag"] - f[0].header["gmag"])
	gr_color.append(f[0].header["gmag"] - f[0].header["rmag"])
	sav.append(f[0].header["v"] - f[0].header["pv"])
	smet.append(f[0].header["Z"])
	smass.append(f[0].header["mass"])

ug_color = np.array(ug_color)
gr_color = np.array(gr_color)
sav      = np.array(sav)
smet     = np.array(smet)
smass    = np.array(smass)

m1 = smass == 0.035507
m2 = smass == 0.34032
m3 = smass == 0.99194
m4 = smass == 3.0883

omass, rmass, otm, rtm, otf, rtf, oz, rz, oav, rav = np.loadtxt("../outputs/photometric_fit/remote_set1/phot_fit.physical", unpack = True)

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-g)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = (smass + 1) * 100, c = sav, color = "w", vmin = 0.0, vmax = 2.0, cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = (smass[m1] + 1) * 100, color = "k")
s2 = ax1.scatter(-99, -99, s = (smass[m2] + 1) * 100, color = "k")
s3 = ax1.scatter(-99, -99, s = (smass[m3] + 1) * 100, color = "k")
s4 = ax1.scatter(-99, -99, s = (smass[m4] + 1) * 100, color = "k")
cb = plt.colorbar(sc, cax = axins1)
cb.set_label("Av[mag]")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2, s3, s4], [0.05, 0.35, 1.00, 3.00], loc = 2, frameon = False, scatterpoints = 1, title = "mass (a/u)")
for i in xrange(4) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("lib.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------

mass_err = err(omass, rmass)
error, sigma = [], []
for i in xrange(134) :
	error.append(np.median(mass_err[i * 100:(i + 1) * 100]))
	sigma.append(percentile(mass_err[i * 100:(i + 1) * 100], 84.0) - percentile(mass_err[i * 100:(i + 1) * 100], 16.0))

error = np.array(error)
sigma = np.array(sigma)

m1 = sigma == sigma.min()
m2 = sigma == sigma.max()

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-r)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = sigma, c = error, vmin = -100, vmax = +100, color = "w", cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = sigma[m1], color = "k")
s2 = ax1.scatter(-99, -99, s = sigma[m2], color = "k")

cb = plt.colorbar(sc, cax = axins1)
cb.set_label("relative error in mass")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2], [round(sigma.min(), 0), round(sigma.max(), 0)], loc = 2, frameon = False, scatterpoints = 1, title = "uncertainty (%)")
for i in xrange(2) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("mass.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------

tm_err = (10 ** (otm - rtm) - 1) * 100.
error, sigma = [], []
for i in xrange(134) :
	error.append(np.median(tm_err[i * 100:(i + 1) * 100]))
	sigma.append(percentile(tm_err[i * 100:(i + 1) * 100], 84.0) - percentile(tm_err[i * 100:(i + 1) * 100], 16.0))

error = np.array(error)
sigma = np.array(sigma)

m1 = sigma == sigma.min()
m2 = sigma == sigma.max()

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-r)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = sigma, c = error, vmin = -100, vmax = +100, color = "w", cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = sigma[m1], color = "k")
s2 = ax1.scatter(-99, -99, s = sigma[m2], color = "k")

cb = plt.colorbar(sc, cax = axins1)
cb.set_label("relative error in MWA")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2], [round(sigma.min(), 0), round(sigma.max(), 0)], loc = 2, frameon = False, scatterpoints = 1, title = "uncertainty (%)")
for i in xrange(2) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("tm.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------

tf_err = (10 ** (otf - rtf) - 1) * 100.
error, sigma = [], []
for i in xrange(134) :
	error.append(np.median(tf_err[i * 100:(i + 1) * 100]))
	sigma.append(percentile(tf_err[i * 100:(i + 1) * 100], 84.0) - percentile(tf_err[i * 100:(i + 1) * 100], 16.0))

error = np.array(error)
sigma = np.array(sigma)

m1 = sigma == sigma.min()
m2 = sigma == sigma.max()

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-r)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = sigma, c = error, vmin = -100, vmax = +100, color = "w", cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = sigma[m1], color = "k")
s2 = ax1.scatter(-99, -99, s = sigma[m2], color = "k")

cb = plt.colorbar(sc, cax = axins1)
cb.set_label("relative error in LWA")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2], [round(sigma.min(), 0), round(sigma.max(), 0)], loc = 2, frameon = False, scatterpoints = 1, title = "uncertainty (%)")
for i in xrange(2) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("tf.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------

av_err = err(oav, rav)
error, sigma = [], []
for i in xrange(134) :
	error.append(np.median(av_err[i * 100:(i + 1) * 100]))
	sigma.append(percentile(av_err[i * 100:(i + 1) * 100], 84.0) - percentile(av_err[i * 100:(i + 1) * 100], 16.0))

error = np.array(error)
sigma = np.array(sigma)

m1 = sigma == sigma.min()
m2 = sigma == sigma.max()

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-r)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = sigma, c = error, vmin = -100, vmax = +100, color = "w", cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = sigma[m1], color = "k")
s2 = ax1.scatter(-99, -99, s = sigma[m2], color = "k")

cb = plt.colorbar(sc, cax = axins1)
cb.set_label("relative error in extinction")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2], [round(sigma.min(), 0), round(sigma.max(), 0)], loc = 2, frameon = False, scatterpoints = 1, title = "uncertainty (%)")
for i in xrange(2) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("av.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------

z_err = err(oz, rz)
error, sigma = [], []
for i in xrange(134) :
	error.append(np.median(z_err[i * 100:(i + 1) * 100]))
	sigma.append(percentile(z_err[i * 100:(i + 1) * 100], 84.0) - percentile(z_err[i * 100:(i + 1) * 100], 16.0))

error = np.array(error)
sigma = np.array(sigma)

m1 = sigma == sigma.min()
m2 = sigma == sigma.max()

plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-r)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

sc = ax1.scatter(ug_color, gr_color, s = sigma, c = error, vmin = -100, vmax = +100, color = "w", cmap = cm.spectral)
s1 = ax1.scatter(-99, -99, s = sigma[m1], color = "k")
s2 = ax1.scatter(-99, -99, s = sigma[m2], color = "k")

cb = plt.colorbar(sc, cax = axins1)
cb.set_label("relative error in metallicity")
xe = np.arange((u - g).min(), (u - g).max(), nbins(u - g)[ - 1])
ye = np.arange((g - r).min(), (g - r).max(), nbins(g - r)[ - 1])
H, xe, ye = np.histogram2d(u - g, g - r, bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], colors = ["#1A1A1A"], linewidths = 2)

leg = ax1.legend([s1, s2], [round(sigma.min(), 0), round(sigma.max(), 0)], loc = 2, frameon = False, scatterpoints = 1, title = "uncertainty (%)")
for i in xrange(2) :
	leg.texts[i].set_fontsize(12)
	leg.legendHandles[i].set_facecolor("none")

axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.savefig("z.png", bbox_inches = "tight")
