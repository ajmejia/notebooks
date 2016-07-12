#!/usr/bin/python

#number of fits read      = 1640320
#number of galaxies found = 906736
#sample size              = 894205

import numpy as np
from matplotlib import cm, rc
import matplotlib.pyplot as plt
from scipy.stats.mstats import scoreatpercentile, mode
from scipy.ndimage.filters import gaussian_filter
import calc_kcor as kcor

rc("legend", numpoints = 1)

def flux(mag) : return 10 ** (- 0.4 * mag)

def binner(x, y, w_sta, nbins, rang = None, ebar = False, per = None) :
	from numpy import array, digitize, lexsort, linspace
	from numpy.ma import average, median

	ind    = lexsort((y, x))
	xs, ys = x[ind], y[ind]

	if rang is None : mn, mx = min(xs), max(xs)
	else            : mn, mx = rang
	
	bins  = linspace(mn, mx, nbins + 1)
	x_cen = (bins[: - 1] + bins[1:])*0.5
	bins  = linspace(mn, mx, nbins)
	ibins = digitize(xs, bins)

	if w_sta   == "median" : y_sta = array([median(ys[ibins == i]) for i in range(1, bins.size + 1)])
	elif w_sta == "mean"   : y_sta = array([average(ys[ibins == i]) for i in range(1, bins.size + 1)])
	elif w_sta == "mode"   : y_sta = array([mode(ys[ibins == i])[0] for i in range(1, bins.size + 1)])

	if ebar   == False                : return x_cen, y_sta
	elif ebar == True and per == None :
		myer = abs(array([scoreatpercentile(ys[ibins == i], 15.8) for i in range(1, bins.size + 1)]) - y_sta)
		pyer = abs(array([scoreatpercentile(ys[ibins == i], 84.0) for i in range(1, bins.size + 1)]) - y_sta)
		yer  = array([myer, pyer])
		return x_cen, y_sta, yer

	elif ebar == True and per != None :
		myer = abs(array([scoreatpercentile(ys[ibins == i], per[0]) for i in range(1, bins.size + 1)]) - y_sta)
		pyer = abs(array([scoreatpercentile(ys[ibins == i], per[1]) for i in range(1, bins.size + 1)]) - y_sta)
		yer = array([myer, pyer])
		return x_cen, y_sta, yer

snr1_sq, snr2_sq, z_spec = np.loadtxt("sample.dat", unpack = True, usecols = (3, 6, 10))

mags   = np.loadtxt("../strateva/data_color_chen_lib_bc03.txt", usecols = (2,3,4,5,6))
colors = np.array([mags[:, j] - mags[:, j + 1] for j in range(np.size(mags, 1) - 1)]).T

#data = np.loadtxt("photo_sample.csv", delimiter = ",", skiprows = 1)
#fiber_u = data[:, 5]
#fiber_g = data[:, 6]
#fiber_r = data[:, 7]
#fiber_i = data[:, 8]
#fiber_z = data[:, 9]
#petro_u = data[:, 25]
#petro_g = data[:, 26]
#petro_r = data[:, 27]
#petro_i = data[:, 28]
#petro_z = data[:, 29]
#z       = data[:, 35]
data = np.loadtxt("output.csv", delimiter = ",", skiprows = 1)
fiber_u = data[:, 2]
fiber_g = data[:, 3]
fiber_r = data[:, 4]
fiber_i = data[:, 5]
fiber_z = data[:, 6]
petro_u = data[:, 12]
petro_g = data[:, 13]
petro_r = data[:, 14]
petro_i = data[:, 15]
petro_z = data[:, 16]
z       = data[:, - 4]
del mags

kcor_fu = kcor.calc_kcor("u", z, "u - r", fiber_u - fiber_r)
kcor_fg = kcor.calc_kcor("g", z, "g - r", fiber_g - fiber_r)
kcor_fr = kcor.calc_kcor("r", z, "g - r", fiber_g - fiber_r)
kcor_fi = kcor.calc_kcor("i", z, "g - i", fiber_g - fiber_i)
kcor_fz = kcor.calc_kcor("z", z, "r - z", fiber_r - fiber_z)
kcor_pu = kcor.calc_kcor("u", z, "u - r", petro_u - petro_r)
kcor_pg = kcor.calc_kcor("g", z, "g - r", petro_g - petro_r)
kcor_pr = kcor.calc_kcor("r", z, "g - r", petro_g - petro_r)
kcor_pi = kcor.calc_kcor("i", z, "g - i", petro_g - petro_i)
kcor_pz = kcor.calc_kcor("z", z, "r - z", petro_r - petro_z)

fiber_u += kcor_fu
fiber_g += kcor_fg
fiber_r += kcor_fr
fiber_i += kcor_fi
fiber_z += kcor_fz
petro_u += kcor_pu
petro_g += kcor_pg
petro_r += kcor_pr
petro_i += kcor_pi
petro_z += kcor_pz
del kcor_fu, kcor_fg, kcor_fr, kcor_fi, kcor_fz, kcor_pu, kcor_pg, kcor_pr, kcor_pi, kcor_pz

#percent = flux(fiber_r) / flux(petro_r) * 100.0
percent = 3.0 / data[:, 24] * 100.0
u_g = fiber_u - fiber_g
g_r = fiber_g - fiber_r
r_i = fiber_r - fiber_i
i_z = fiber_i - fiber_z

hs = (0.02 <= z) & (z <= 0.03)
z1 = (0.1 <= z) & (z <= 0.2)
z2 = (0.2 <= z) & (z <= 0.3)
z3 = (0.3 <= z) & (z <= 0.4)
zx = ~hs & ~z1 & ~z2 & ~z3

col = ["gray", "#FFA500", "#1E90FF", "#CB7E9C", "#9ABD5C"]

fig = plt.figure(figsize = (8, 8))
fig.subplots_adjust(left = 0.09, bottom = 0.07, right = 0.96, top = 0.97, hspace = 0.11, wspace = 0.11)

ax1 = plt.subplot2grid((4, 4), (0, 0), rowspan = 3, colspan = 3)
ax2 = plt.subplot2grid((4, 4), (3, 0), colspan = 3)
ax3 = plt.subplot2grid((4, 4), (0, 3), rowspan = 3)

ax1.plot(z[zx], percent[zx], "s", mfc = col[0], mec = "none", ms = 1)
ax1.plot(z[hs], percent[hs], "s", mfc = col[1], mec = "none", ms = 1, label = "0.02 - 0.03")
ax1.plot(z[z1], percent[z1], "s", mfc = col[2], mec = "none", ms = 1, label = "0.10 - 0.20")
ax1.plot(z[z2], percent[z2], "s", mfc = col[3], mec = "none", ms = 1, label = "0.20 - 0.30")
ax1.plot(z[z3], percent[z3], "s", mfc = col[4], mec = "none", ms = 1, label = "0.30 - 0.40")
x, y, yer = binner(z, percent, "median", 50, (0.01, 0.5), True)
ax1.plot(x, y, "--k", lw = 1.5)

leg = ax1.legend(loc = 4, markerscale = 10, fontsize = 12, frameon = False, title = "redshift cuts")
for i, line in enumerate(leg.get_lines()) : line.set_alpha(1.0)

ax1.tick_params(labelbottom = False)
ax1.set_xlim(0.0, 0.6)
ax1.set_ylim(0.0, 100.0)
ax1.set_ylabel("flux percent in r band")

density, zi, patches = ax2.hist(z, 50, ec = col[0], fc = "none", normed = True, range = (0.0, 0.6), hatch = "////", lw = 1.5)
ax2.set_xlabel("redshift")

bwidth = zi[1] - zi[0]

for patch in patches :
	if 0.02 <= patch.get_x() + bwidth * 0.5 <= 0.03 : patch.set_ec(col[1])
	elif 0.1 <= patch.get_x() + bwidth * 0.5 <= 0.2 : patch.set_ec(col[2])
	elif 0.2 <= patch.get_x() + bwidth * 0.5 <= 0.3 : patch.set_ec(col[3])
	elif 0.3 <= patch.get_x() + bwidth * 0.5 <= 0.4 : patch.set_ec(col[4])

data_list = [percent[hs], percent[z1], percent[z2], percent[z3]]
ax3.hist(data_list, 50, color = col[1:], ec = "none", normed = True, orientation = "horizontal", range = (0.0, 100.0), rwidth = 1.0)

ax3.tick_params(labelleft = False)
ax3.set_xticks(np.arange(0.0, 0.07, 0.03))

plt.savefig("flux_percent_z", bbox_inches = "tight")

fig = plt.figure(figsize = (6, 8))
fig.subplots_adjust(left = 0.12, bottom = 0.07, right = 0.91, top = 0.97, wspace = 1.35, hspace = - 5.39)

ax1 = plt.subplot2grid((16, 16), (0, 0), rowspan = 8, colspan = 15)
ax2 = plt.subplot2grid((16, 16), (8, 0), rowspan = 8, colspan = 15)
ax3 = plt.subplot2grid((16, 16), (0, 15), rowspan = 8)
ax4 = plt.subplot2grid((16, 16), (8, 15), rowspan = 8)

xe = np.arange(colors[:, 0].min(), colors[:, 0].max(), 0.02)
ye = np.arange(colors[:, 1].min(), colors[:, 1].max(), 0.01)
H, xe, ye = np.histogram2d(colors[:, 0], colors[:, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
cont = ax1.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r)
cb = plt.colorbar(cont, cax = ax3)
cb.set_label("models per bin")

H, xe, ye = np.histogram2d(u_g[hs], g_r[hs], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [1], colors = col[1], linewidths = 2)

H, xe, ye = np.histogram2d(u_g[z1], g_r[z1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [1], colors = col[2], linewidths = 2)

H, xe, ye = np.histogram2d(u_g[z2], g_r[z2], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [1], colors = col[3], linewidths = 2)

H, xe, ye = np.histogram2d(u_g[z3], g_r[z3], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [1], colors = col[4], linewidths = 2)

ax1.set_xlim(0.8, 2.4)
ax1.set_ylim(0.2, 1.7)
ax1.set_xlabel("(u-g)[mag]")
ax1.set_ylabel("(g-r)[mag]")

xe = np.arange(colors[:, 2].min(), colors[:, 2].max(), 0.02)
ye = np.arange(colors[:, 3].min(), colors[:, 3].max(), 0.01)
H, xe, ye = np.histogram2d(colors[:, 2], colors[:, 3], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
cont = ax2.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r)
cb = plt.colorbar(cont, cax = ax4)
cb.set_label("models per bin")

H, xe, ye = np.histogram2d(r_i[hs], i_z[hs], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax2.contour(xe[:- 1], ye[:- 1], H.T, levels = [4], colors = col[1], linewidths = 2)

H, xe, ye = np.histogram2d(r_i[z1], i_z[z1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax2.contour(xe[:- 1], ye[:- 1], H.T, levels = [4], colors = col[2], linewidths = 2)

H, xe, ye = np.histogram2d(r_i[z2], i_z[z2], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax2.contour(xe[:- 1], ye[:- 1], H.T, levels = [4], colors = col[3], linewidths = 2)

H, xe, ye = np.histogram2d(r_i[z3], i_z[z3], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax2.contour(xe[:- 1], ye[:- 1], H.T, levels = [4], colors = col[4], linewidths = 2)

ax2.set_xlim(0.05, 0.75)
ax2.set_ylim(0.0, 0.6)
ax2.set_xlabel("(r-i)[mag]")
ax2.set_ylabel("(i-z)[mag]")

plt.savefig("color_sampling", bbox_inches = "tight")
