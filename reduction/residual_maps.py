#!/usr/bin/python

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.stats as st
import pyfits as pyf

from matplotlib.mlab import griddata
from scipy.ndimage.filters import gaussian_filter
from thesis_redux_tools import nbins, binned_stat, err, binner

# READ DATA ========================================================================================

# read sample properties ---------------------------------------------------------------------------

sfh_name_sam = np.loadtxt("../../inputs/SFHs_set1/set1.log", dtype = "|S")

u_sam, g_sam, r_sam, i_sam, z_sam       = [], [], [], [], []
#tm_sam, tf_sam, Z_sam, Av_sam, mass_sam = [], [], [], [], []
phy_sam = []
for name in sfh_name_sam :
	f = pyf.open("../../inputs/SFHs_set1/" + name)
	
	u_sam.append(f[0].header["umag"])
	g_sam.append(f[0].header["gmag"])
	r_sam.append(f[0].header["rmag"])
	i_sam.append(f[0].header["imag"])
	z_sam.append(f[0].header["zmag"])
  
	phy_sam.append([f[0].header["mass"], f[0].header["mwla"], f[0].header["rfwla"], f[0].header["z"], f[0].header["v"] - f[0].header["pv"]])

u_sam = np.array(u_sam)
g_sam = np.array(g_sam)
r_sam = np.array(r_sam)
i_sam = np.array(i_sam)
z_sam = np.array(z_sam)

phy_sam = np.array(phy_sam)
phy_mod = np.loadtxt("../../outputs/photometric_fit/remote_set1/phot_fit.models", usecols = range(5, 10))

med = binned_stat(phy_mod, 100, stat = "median")
std = binned_stat(phy_mod, 100, stat = "stdev")

# read SFH library data ----------------------------------------------------------------------------

sfhs_cat = "../../inputs/run07_sfh_catalog.txt"
u_lib, g_lib, r_lib, i_lib, z_lib, V, pV, mu = np.genfromtxt(sfhs_cat, missing = '""', usecols = range(39, 44) + [38, 25, 14], unpack = True)

# define mask --------------------------------------------------------------------------------------

dmask = (V - pV < 2.0) & (mu <= 1.0)

# MAKE PLOT ========================================================================================

# define color-color plane -------------------------------------------------------------------------

lcolor_x = (u_lib - g_lib)[dmask]
lcolor_y = (g_lib - r_lib)[dmask]
#lcolor_y = (r_lib - i_lib)[dmask]

scolor_x = (u_sam - g_sam)
scolor_y = (g_sam - r_sam)
#scolor_y = (r_sam - i_sam)

# define plot properties ---------------------------------------------------------------------------

fig, axs = plt.subplots(2, 2, sharex = "col", sharey = "row", figsize = (10,9))

# construct 0.1 models per bin contour of SSAG -----------------------------------------------------

nx, xi, xf, bsx = nbins(lcolor_x)
ny, yi, yf, bsy = nbins(lcolor_y)

print("bin sizes:")
print("xbins : {0:5.3f}".format(bsx))
print("ybins : {0:5.3f}".format(bsy))

H, xe, ye = np.histogram2d(lcolor_x, lcolor_y, bins = (nx, ny), normed = True)
H         = gaussian_filter(H, sigma = 2)

# plot contour and residual of sample model galaxies -----------------------------------------------

my_cm = mpl.colors.ListedColormap(np.loadtxt("colormap.dat") / 255)
for i, ax in enumerate(np.ravel(axs)) :
	ax.set_xlim(0.5, 2.5)
	ax.set_ylim(0.0, 1.25)

	ax.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.1], colors = ["#1A1A1A"], lw = 2.5)

	sizes  = std[:, i] * 100
	if i == 0 or i == 3 : colors = err(phy_sam[:, i], med[:, i])
	else                : colors = err(phy_sam[:, i], med[:, i], relative = False)

	xg = np.linspace(0.5, 2.5, 100)
	yg = np.linspace(0.0, 1.25, 100)
	zg = griddata(scolor_x, scolor_y, std[:, i], xg, yg, interp = "linear")

	#cont = ax.contourf(xg, yg, zg, vmin = -1.0, vmax = +1.0, levels = np.linspace(-1., +1., 9), cmap = my_cm)
	cont = ax.contour(xg, yg, zg, 3)
	plt.colorbar(cont)

	#if i == 0 :
		#sc = ax.scatter(scolor_x, scolor_y, lw = 0, c = colors, s = sizes, vmin = -1.0, vmax = +1.0, cmap = my_cm)
		#cb = plt.colorbar(sc)
		#cb.set_label("residual")
		#cb.set_ticks(np.linspace(-1., +1., 9))
	#else : ax.scatter(scolor_x, scolor_y, lw = 0, c = colors, s = sizes, vmin = -1.0, vmax = +1.0, cmap = my_cm)

plt.tight_layout()
plt.savefig("../../plots/photometric_fit/residual_maps.png", bbox_inches = "tight")
plt.show()
