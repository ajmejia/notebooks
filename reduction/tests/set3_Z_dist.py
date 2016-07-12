from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np

rc("text", usetex = True)

Z = np.loadtxt("../../inputs/set3_catalog.txt", usecols = (12,))

plt.hist(Z, 10, histtype = "step", ec = "#7F7F7F", lw = 2, hatch = "/")

for z in [0.004, 0.008, 0.02, 0.05]: plt.axvline(z/0.02, ls = "--", lw = 2, color = "#1E90FF")

plt.xlabel("$Z/Z_\odot$", size = 16)
plt.ylabel("Counts", size = 16)

plt.xlim(0.0, 2.6)
plt.ylim(0, 17)

plt.tight_layout()

plt.savefig("../../plots/tests/set3_Z2_tex.eps")
plt.show()
