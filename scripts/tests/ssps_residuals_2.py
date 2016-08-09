#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import sys


if len(sys.argv) > 1 : fid = sys.argv[1]
else                 : fid = "default"

t_i, t_j, age_i, age_j, theta = np.loadtxt("../../outputs/tests/ssps_comp_" + fid + ".txt", unpack = True)

ax = plt.subplot(111, aspect = "equal")
sc = ax.scatter(t_i, t_j, lw = 0, s = 7, c = theta, vmin = 0, vmax = 90)
cb = plt.colorbar(sc)
lm = min(t_i), max(t_i)
fs = 16
tk = range(1, 196, 20)
tl = map(lambda fig : "{0:15.2e}".format(fig), age_j[:195][tk])

ax.set_xlim(lm)
ax.set_ylim(lm)
cb.set_label(r"$\theta$", fontsize = fs)
ax.set_title(fid.replace("_", " "))
ax.set_xlabel("age", fontsize = fs)
ax.set_ylabel("age", fontsize = fs)

ax.set_xticks(tk)
ax.set_yticks(tk)
ax.set_xticklabels(tl, rotation = 45)
ax.set_yticklabels(tl)

plt.subplots_adjust(right = 0.75, bottom = 0.15)
plt.savefig("../../plots/ssps_com_" + fid, bbox_inches = "tight")
plt.show()
