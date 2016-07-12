import numpy as np
import matplotlib.pyplot as plt

SDSS   = np.loadtxt("../../inputs/SDSS_filters.txt")
COSMOS = np.loadtxt("../../inputs/COSMOS_filters.txt")
JPAS   = np.loadtxt("../../inputs/J-PAS_filters.txt")

plt.semilogx(SDSS[:, 0], SDSS[:, 1], "--k", lw = 2, label = "SDSS")
plt.semilogx(COSMOS[:, 0], COSMOS[:, 1], "-", color = "#7F7F7F", lw = 1, label = "COSMOS")
plt.semilogx(JPAS[:, 0], JPAS[:, 1], "-", color = "#FF4500", lw = 0.5, label = "J-PAS")
plt.legend(loc = 2, frameon = False, numpoints = 1)
plt.xlabel("wavelength[A]")
plt.ylabel("transmision[%]")
#plt.xlim(2000, 51000)
plt.ylim(0.0001, 1)
plt.tight_layout()
plt.show()
