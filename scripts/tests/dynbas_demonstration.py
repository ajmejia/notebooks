from pylab import *

wl, fl, fl3d = loadtxt("dynbasfit_processed_fullsed_spSpec-52146-0658-247.log", usecols = (0, 1, 5), unpack = True)
wlm, g1, g2  = loadtxt("gen_1_2", unpack = True)
wlm, g3      = loadtxt("gen_3", unpack = True)

coeff = array([44082384.0, 1.22554867E+09, 3065331.50])

#figure(figsize = (12, 6))

ax1 = subplot2grid((3, 3), (0, 0), rowspan = 3, colspan = 2)
ax2 = subplot2grid((3, 3), (0, 2))
ax3 = subplot2grid((3, 3), (1, 2))
ax4 = subplot2grid((3, 3), (2, 2))

ax1.set_xlabel(r"$\lambda$ (\AA)", fontsize = 15)
ax1.set_ylabel(r"$\text{L}_\lambda$ ($10^5\,$L$\odot$/\AA)", fontsize = 15)

ax4.set_xlabel(r"$\lambda$ (\AA)", fontsize = 15)

for ax in [ax2, ax3] : ax.set_xticklabels([])
for ax in [ax2, ax3, ax4] : ax.set_xticks([3000, 5000, 7000, 9000])

ax1.plot(wl, fl   / 1e5, "-k", lw = 1, label = "target")
ax1.plot(wl, fl3d / 1e5, "-", color = "#80CC28", lw = 1, label = "model")
ax1.legend(loc = 1, prop = {"size" : 15})

ax2.plot(wlm, g1 * coeff[0] / 1e5, "-", color = "#266BBD", lw = 1)
ax3.plot(wlm, g2 * coeff[1] / 1e5, "-", color = "#ED0E58", lw = 1)
ax4.plot(wlm, g3 * coeff[2] / 1e5, "-", color = "#F68512", lw = 1)

ax2.text(5000, 1.3, r"$M_\star \%={0:.2f}$".format(coeff[0] / coeff.sum() * 100.0), fontsize = 15)
ax3.text(5000, 0.3, r"$M_\star \%={0:.2f}$".format(coeff[1] / coeff.sum() * 100.0), fontsize = 15)
ax4.text(5000, 0.05, r"$M_\star \%={0:.2f}$".format(coeff[2] / coeff.sum() * 100.0), fontsize = 15)

tight_layout()
savefig("dynbas_demonstration.pdf", bbox_inches = "tight")
show()
