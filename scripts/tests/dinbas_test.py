#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

masses_new, masses_old = np.loadtxt("fort.1")
m2dn = masses_new[:2]
m3dn = masses_new[2:]
m2do = masses_old[:2]
m3do = masses_old[2:]

wl, fl, f2d1, f2d2, f3d1, f3d2, f3d3 = np.loadtxt("fort.2")
fit_new2d = m2dn[0] * f2d1 + m2dn[1] * f2d2
fit_new3d = m3dn[0] * f3d1 + m3dn[1] * f3d2 + m3dn[2] * f3d3

wl, fl, f2d1, f2d2, f3d1, f3d2, f3d3 = np.loadtxt("fort.3")
fit_old2d = m2do[0] * f2d1 + m2do[1] * f2d2
fit_old3d = m3do[0] * f3d1 + m3do[1] * f3d2 + m3do[2] * f3d3


ax1 = plt.subplot(211)
ax2 = plt.subplot(212)

ax1.plot(wl, fl, color = "k")
ax1.plot(wl, fit_old3d, color = "#FFA500")

ax2.plot(wl, fl, color = "k")
ax2.plot(wl, fit_new3d, color = "#FFA500")

plt.show()
