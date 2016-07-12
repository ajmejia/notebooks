import numpy as np
from matplotlib import rc, cm
import matplotlib.pyplot as plt
import pyfits as pyf
from scipy.ndimage.filters import gaussian_filter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

data   = np.loadtxt("../../TEG/data/photometry/total_photometry.txt", usecols = range(1,4) + [11])
mags   = data[:, : - 1]
z      = data[:, -1]
colors = np.array([mags[:, j] - mags[:, j + 1] for j in xrange(np.size(mags, 1) - 1)]).T

table = np.genfromtxt("../inputs/run07_sfh_catalog.txt", missing = '""', usecols = range(39, 42) + [14, 25, 38])
mu = table[:, - 3]
Av = table[:, - 1] - table[:, - 2]
mags = table[:, :3]
colors_m = np.array([mags[:, j] - mags[:, j + 1] for j in range(np.size(mags, 1) - 1)]).T

del table, mags, data

fitsnames = np.loadtxt("../inputs/SFHs_set1/set1.log", dtype = "|S")

ug_color, gr_color = [], []
mu_s = []
for name in fitsnames :
	
	f = pyf.open("../inputs/SFHs/" + name)
	
	ug_color.append(f[0].header["umag"] - f[0].header["gmag"])
	gr_color.append(f[0].header["gmag"] - f[0].header["rmag"])
	mu_s.append(f[0].header["mu"])

fig = plt.figure(figsize = (7, 6))
ax1 = plt.subplot(111, xlim = (0.5, 2.5), ylim = (0.0,1.25), xlabel = "(u-g)[mag]", ylabel = "(g-r)[mag]")
axins1 = inset_axes(ax1, width="2.5%", height="40%", loc=4)

mask  = z < 0.04
xe = np.arange(colors[mask, 0].min(), colors[mask, 0].max(), 0.02)
ye = np.arange(colors[mask, 1].min(), colors[mask, 1].max(), 0.01)
H, xe, ye = np.histogram2d(colors[mask, 0], colors[mask, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
cont = ax1.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r, levels = np.arange(0.0, 7.0, 1.0))
cb = plt.colorbar(cont, cax = axins1)
cb.set_label("galaxies per bin")

#mask = Av < 2.0
mask = (mu < 1.0) & (Av < 2.0)
H, xe, ye = np.histogram2d(colors_m[mask, 0], colors_m[mask, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
ax1.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], linewidths = 2.5, colors = ["#1E90FF"])

mu_s = np.array(mu_s)
ug_color = np.array(ug_color)[mu_s<=1.0]
gr_color = np.array(gr_color)[mu_s<=1.0]
s, = ax1.plot(ug_color, gr_color, "x", mfc = "none", mec = "orangered", mew = 1.2)
ax1.plot([0.5, 2.5], [-0.5+2.22, -2.5+2.22], "--", lw = 1.5, color = "#1A1A1A")
#r, = ax1.plot(ug_color[:16], gr_color[:16], "x", mfc = "none", mec = "gold", mew = 1.2)


axins1.tick_params(right = False, left = True, direction = "out", labelright = "off", labelleft = "on")
axins1.yaxis.set_label_position("left")

plt.draw()
#plt.savefig("localrun_lib.png", bbox_inches = "tight")
plt.show()
