#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.mlab import griddata
from matplotlib import cm, rc

rc("text", usetex=True)

if len(sys.argv) > 1 : fid = sys.argv[1]
else                 : fid = "default"

age_i, age_j, theta = np.loadtxt("../../outputs/tests/ssps_comp_" + fid + ".txt", unpack=True, usecols=(2, 3, 4))

theta[np.isnan(theta)] = 0.0

lim = (age_i.min(), age_i.max())
t_i = np.logspace(np.log10(lim[0]), np.log10(lim[1]), 50)
t_j = np.logspace(np.log10(lim[0]), np.log10(lim[1]), 50)
theta_ij = griddata(age_i, age_j, theta, t_i, t_j)

ax = plt.subplot(111, aspect="equal", xscale="log", yscale="log")
ct = ax.contourf(t_i, t_j, theta_ij, levels=range(0, 100, 10), cmap=cm.rainbow)
cb = plt.colorbar(ct)

ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_title(fid.replace("_", " "), fontsize=16)
ax.set_xlabel(r"$t_i$", fontsize=14)
ax.set_ylabel(r"$t_j$", fontsize=14)
cb.set_label(r"$\theta_{ij}$", fontsize=14)

plt.savefig("../../plots/ssps_com_" + fid, bbox_inches = "tight")
plt.show()
