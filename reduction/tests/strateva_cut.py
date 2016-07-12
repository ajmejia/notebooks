from matplotlib import rc, cm, use
#use("Qt4Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage.filters import gaussian_filter
from thesis_redux_tools import *

#rc('font', **{'family':'serif', 'serif':['Roman']})
rc('text', usetex = True)

# OLD LIBRARY ======================================================================================
data   = read_mybin("../../../TEG/output/thesis_data/rest_frame_run/library/for_dinbas.bin", to_return = ["fluxes"], kind = "empi")
mags   = - 2.5*np.log10(array(data["fluxes"])[0]) - 48.60
colors_m_old = np.array([mags[:, j] - mags[:, j + 1] for j in range(np.size(mags, 1) - 1)]).T
# ==================================================================================================

# NEW LIBRARY ======================================================================================
#mags = np.loadtxt("data_color_chen_lib.txt", usecols = (8,9,10,11,12))
#mags = np.loadtxt("data_color_chen_lib.txt", usecols = (3,4,5,6,7))
#mags = np.loadtxt("data_color_chen_lib_bc03.txt", usecols = (7,8,9,10,11))
#mags = np.loadtxt("data_color_chen_lib_bc03.txt", usecols = (2,3,4,5,6))
#mags = np.loadtxt("table.txt", usecols = (7,8,9,10,11)) # plain
#mags = np.loadtxt("table.txt", usecols = (2,3,4,5,6)) # dusty
table = np.genfromtxt("../../inputs/run07_sfh_catalog.txt", missing = '""', usecols = range(39, 44) + range(26, 31) + [25, 38])
Av = table[:, - 1] - table[:, - 2]
mags = table[:, :5]
colors_m_new = np.array([mags[:, j] - mags[:, j + 1] for j in range(np.size(mags, 1) - 1)]).T
# ==================================================================================================

data = np.loadtxt("../../../TEG/data/photometry/total_photometry.txt", usecols = range(1,6) + [11])
mags = data[:, : - 1]
z    = data[:, -1]
mask = (z < 0.04)
colors = np.array([mags[mask, j] - mags[mask, j + 1] for j in xrange(np.size(mags, 1) - 1)]).T
del data, mags

labels = [r"($u-g$)[mag]", r"($g-r$)[mag]", r"($r-i$)[mag]", r"($i-z$)[mag]"]

def strateva_cut(x): return - x + 2.22

fig = plt.figure(figsize = (7, 6))
#plt.subplot(111, xlim = (0.33, 2.35), ylim = (-0.074,1.25))
xe = np.arange(colors[:, 0].min(), colors[:, 0].max(), 0.02)
ye = np.arange(colors[:, 1].min(), colors[:, 1].max(), 0.01)
H, xe, ye = np.histogram2d(colors[:, 0], colors[:, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
plt.contourf(xe[:- 1], ye[:- 1], H.T, cmap = cm.gray_r, levels = np.arange(0.0, 7.0, 1.0))
cb = plt.colorbar()

H, xe, ye = np.histogram2d(colors_m_old[:, 0], colors_m_old[:, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
plt.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], linewidths = 2.5, colors = ["#FF6347"])

H, xe, ye = np.histogram2d(colors_m_new[:, 0], colors_m_new[:, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
plt.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], linewidths = 2.5, colors = ["#5F9EA0"])

mask = Av < 1.5
H, xe, ye = np.histogram2d(colors_m_new[mask, 0], colors_m_new[mask, 1], bins = (xe, ye), normed = True)
H = gaussian_filter(H, sigma = 2)
plt.contour(xe[:- 1], ye[:- 1], H.T, levels = [0.2], linewidths = 2.5, colors = ["#FFA500"])

plt.xlabel(labels[0], fontsize = 16)
plt.ylabel(labels[1], fontsize = 16)
cb.set_label(r"counts per bin", fontsize = 16)

#plt.savefig("strateva_sep_old.pdf", bbox_inches = "tight")

plt.show()
