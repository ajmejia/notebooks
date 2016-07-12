import numpy as np
from matplotlib import use
use("Qt4Agg")
import matplotlib.pyplot as plt
from matplotlib import rc

#rc("text", usetex = True)

omass, rmass, otm, rtm, otf, rtf, oz, rz, oav, rav = np.loadtxt("../outputs/phot_fit.physical", unpack = True)

plt.figure()
plt.subplot(111, xlim = (0, 5), ylim = (0, 5), xlabel = r"$M_\star$ simulada", ylabel = r"$M_\star$ recuperada", aspect = "equal")
plt.plot(omass, rmass, ",", color = "gray")
plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)
plt.savefig("photo_mass_ret.png", bbox_inches = "tight")

plt.figure()
plt.subplot(111, xlim = (0, 13.75e9), ylim = (0, 13.75e9), xlabel = r"$\left<t\right>_M$ simulada", ylabel = r"$\left<t\right>_M$ recuperada")
plt.plot(otm, rtm, ",", color = "gray")
plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)

plt.figure()
plt.subplot(111, xlim = (0, 13.75e9), ylim = (0, 13.75e9), xlabel = r"$\left<t\right>_L$ simulada", ylabel = r"$\left<t\right>_L$ recuperada")
plt.plot(otf, rtf, ",", color = "gray")
plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)

plt.figure()
plt.subplot(111, xlim = (0, 6), ylim = (0, 6), xlabel = r"$A_V$ simulada", ylabel = r"$A_V$ recuperada")
plt.plot(oav, rav, ",", color = "gray")
plt.plot((0, 6), (0, 6), "--", color = "r", lw = 1.5)

plt.figure()
plt.subplot(111, xlim = (0, 6), ylim = (0, 6), xlabel = r"$A_V$ simulada", ylabel = r"$A_V$ recuperada")
plt.plot(oz, rz, ",", color = "gray")
plt.plot((0, 6), (0, 6), "--", color = "r", lw = 1.5)

plt.show()
