#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as st
import data_loader as dl
from itertools import product
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib import cm
import thesis_redux_tools as rt

name_pho  = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", usecols = (0,), dtype = "|S")
name_sed  = np.loadtxt("../../outputs/spectroscopic_fit/names.log", dtype = "|S")
order_sed = [j for i in xrange(120) for j in xrange(12000) if name_pho[i] == name_sed[j][:12] + ".fits.gz"]

osed = np.loadtxt("../../outputs/spectroscopic_fit/table_din.v3.log")[order_sed]

sdss = dl.load_data("SDSS")
jpas = dl.load_data("JPAS")

nams = ["M_lib", "M_mod", "log_t_M_lib", "log_t_M_mod", "log_t_L_lib", "log_t_L_mod", "Z_M_lib", "Z_M_mod", "Av_lib", "Av_mod"]

sdss.physical[nams[6]] = np.log10(sdss.physical[nams[6]])
res_sed = [rt.err(sdss.physical[n], osed[:, i]) if i == 0 else rt.err(sdss.physical[n], osed[:, i], False) for i, n in enumerate(nams[::2])]
sdss.physical[nams[6]] = 10 ** sdss.physical[nams[6]]

fig, axs = plt.subplots(2, 3, sharex = True, figsize = (10, 6))
axins    = inset_axes(axs[0, 0], width = "50%", height = "4%", loc = 2)

for i, j in product(xrange(axs.shape[0]), xrange(axs.shape[1])) :

  if j == 0 :
    if i == 0 :
      axs[i, j].set_title("SDSS")
      axs[i, j].set_ylim(-1.5, 3.5)
      axs[i, j].set_ylabel(r"$\Delta\,M/M_\textrm{SSAG}$", fontsize = 15)
      sc = axs[i, j].scatter(sdss.residuals["log_t_M"], sdss.residuals["M"], lw = 0, s = 7, c = sdss.residuals["Av"], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(sdss.residuals["log_t_M"], sdss.residuals["M"], "mean", 10)
      axs[i, j].plot(x, y, "--k")
    if i == 1 :
      axs[i, j].set_ylim(-1.5, 1.5)
      axs[i, j].set_xlabel(r"\Delta\,\left<\log{t}\right>_M", fontsize = 15)
      axs[i, j].set_ylabel(r"$\Delta\,\left<\log{Z/Z_\odot}\right>_M$", fontsize = 15)
      axs[i, j].scatter(sdss.residuals["log_t_M"], sdss.residuals["log_Z_M"], lw = 0, s = 7, c = sdss.residuals["Av"], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(sdss.residuals["log_t_M"], sdss.residuals["log_Z_M"], "mean", 10)
      axs[i, j].plot(x, y, "--k")

  elif j == 1 :
    axs[i, j].set_yticklabels([])
    axs[i, j].set_xticks(axs[i, j].get_xticks())

    if i == 0 :
      axs[i, j].set_title("J-PAS")
      axs[i, j].set_ylim(-1.5, 3.5)
      axs[i, j].scatter(jpas.residuals["log_t_M"], jpas.residuals["M"], lw = 0, s = 7, c = jpas.residuals["Av"], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(jpas.residuals["log_t_M"], jpas.residuals["M"], "mean", 10)
      axs[i, j].plot(x, y, "--k")
    if i == 1 :
      axs[i, j].set_xlabel(r"$\Delta\,\left<\log{t}\right>_M$", fontsize = 15)
      axs[i, j].set_ylim(-1.5, 1.5)
      axs[i, j].scatter(jpas.residuals["log_t_M"], jpas.residuals["log_Z_M"], lw = 0, s = 7, c = jpas.residuals["Av"], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(jpas.residuals["log_t_M"], jpas.residuals["log_Z_M"], "mean", 10)
      axs[i, j].plot(x, y, "--k")

  elif j == 2 :
    axs[i, j].set_yticklabels([])
    axs[i, j].set_xticks(axs[i, j].get_xticks())

    if i == 0 :
      axs[i, j].set_title("Spectroscopy")
      axs[i, j].set_ylim(-1.5, 3.5)
      axs[i, j].scatter(res_sed[1], res_sed[0], lw = 0, s = 7, c = res_sed[4], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(res_sed[1], res_sed[0], "mean", 10)
      axs[i, j].plot(x, y, "--k")
    if i == 1 :
      axs[i, j].set_xlabel(r"$\Delta\,\left<\log{t}\right>_M$", fontsize = 15)
      axs[i, j].set_ylim(-1.5, 1.5)
      axs[i, j].scatter(res_sed[1], res_sed[3], lw = 0, s = 7, c = res_sed[4], vmin = sdss.residuals["Av"].min(), vmax = sdss.residuals["Av"].max())
      x, y = rt.binner(res_sed[1], res_sed[3], "mean", 10)
      axs[i, j].plot(x, y, "--k")

  axs[i, j].axvline(ls = ":", color = "#FF4500")
  axs[i, j].axhline(ls = ":", color = "#FF4500")


axs[0,0].set_xlim(-2,1.5)
cb = plt.colorbar(sc, cax = axins, orientation = "horizontal")
cb.set_ticks(np.linspace(-2.4, +2.4, 5))
cb.set_label(r"$\Delta\,A_V$", fontsize = 12)

plt.subplots_adjust(hspace = 0.07, wspace = 0.05)
plt.tight_layout()
plt.savefig("/home/alfredo/Dropbox/paper/residuals_coor.png", bbox_inches = "tight")
plt.show()
