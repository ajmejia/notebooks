import numpy as np
from matplotlib import use
use("GTKAgg")
import matplotlib.pyplot as plt

omass, rmass, otm, rtm, otf, rtf, oav, rav = np.loadtxt("spec_fit", unpack = True)

#plt.figure()
#plt.subplot(111, xlim = (0, 4), ylim = (0, 5), xlabel = r"$M_\star$ observada", ylabel = r"$M_\star$ recuperada")
#plt.plot(omass, rmass, "ok")
#plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)
#plt.savefig("spec_mass_ret", bbox_inches = "tight")

#plt.figure()
#plt.subplot(111, xlim = (0, 4), ylim = (0, 5), xlabel = r"$M_\star$ observada", ylabel = r"$M_\star$ recuperada")
#plt.plot(otm, rtm, "ok")
#plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)

#plt.figure()
#plt.subplot(111, xlim = (0, 4), ylim = (0, 5), xlabel = r"$M_\star$ observada", ylabel = r"$M_\star$ recuperada")
#plt.plot(otf, rtf, "ok")
#plt.plot((0, 5), (0, 5), "--", color = "r", lw = 1.5)
#
#plt.figure()
#plt.subplot(111, xlim = (0, 6), ylim = (0, 6), xlabel = r"$A_V$ observada", ylabel = r"$A_V$ recuperada")
#plt.plot(oav, rav, "ok")
#plt.plot((0, 6), (0, 6), "--", color = "r", lw = 1.5)
#plt.savefig("spec_Av_ret", bbox_inches = "tight")

plt.show()
