import numpy as np
import matplotlib.pyplot as plt
import pres_style, os

from PyTools import binner

plt.ion()

ssag = np.genfromtxt(os.path.expandvars("$master/inputs/SSAG_catalog.txt"), dtype=None, names=True)
rgal = np.genfromtxt("data/RealDeal/catalog.txt", dtype=None, names=True)

Av_ssag = ssag["V"]-ssag["pV"]
Av_rgal = np.loadtxt("data/dynbas_output.log", usecols=(8,))

ur_ssag = ssag["umag"]-ssag["rmag"]
ur_rgal = rgal["u_mag"]-rgal["r_mag"]

xm, ym = binner(ur_ssag, Av_ssag, "median", 10, rang=(0.0, 3.5))
xo, yo = binner(ur_rgal, ur_rgal, "median", 10, rang=(0.0, 3.5))

plt.plot(xm, ym, "-", label="SSAG")
plt.plot(xo, yo, "-", label="Real Gals.")
plt.legend(loc=0, fontsize=10)

plt.xlabel("Rest frame u-r (mag)")
plt.ylabel("Internal Extinction (mag)")

#plt.savefig("reddening_ssag_vs_real", bbox_inches="tight")

