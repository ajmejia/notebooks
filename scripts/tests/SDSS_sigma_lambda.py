import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as inter

SDSS   = np.loadtxt("../../inputs/SDSS_filters.txt")
JPAS   = np.loadtxt("../../inputs/J-PAS_filters.txt")
table  = np.loadtxt("../../../TEG/data/photometry/total_photometry.txt", usecols = range(1, 12))
#table  = np.loadtxt("../../inputs/DR8/DR8.csv", usecols = range(17, 17 + 10) + [5], skiprows = 1, delimiter = ",")

mask_z = table[:, -1] < 0.03
mask_m = (table[:, 5] > 0)&(table[:, 6] > 0)&(table[:, 7] > 0)&(table[:, 8] > 0)&(table[:, 9] > 0)
mask_r = (table[:, 2] > 12.)&(table[:, 2] < 18.)&(table[:, 0] > 12.)&(table[:, 0] < 18.)&(table[:, 4] > 12.)&(table[:, 4] < 18.)
mask   = (mask_z)&(mask_m)&(mask_r)
table  = table[mask, :]
wlength= np.array([3.5618E+03, 4.7189E+03, 6.1852E+03, 7.4997E+03, 8.9615E+03])
fluxes = 10 ** (0.4 * (22.5 - table[:, :5])) * 3.631e-6 * 1e-23 * 3e18 / wlength ** 2
sigmas = fluxes * table[:, 5:10] / (2.5 * np.log10(np.exp(1.0)))
#fluxes = table[:, :5]  #* 3.631e-6 * 1e-23 * 3e18 / wlength ** 2
#sigmas = 1/np.sqrt(table[:, 5:10]) #* 3.631e-6 * 1e-23 * 3e18 / wlength ** 2
wl_jpas = np.array([3537.4, 3789.4, 3902.7, 4002.4, 4102.0, 4201.8, 4301.6, 4401.5, 4501.2, 4600.8,
                    4700.8, 4800.7, 4900.7, 5000.4, 5100.2, 5200.2, 5300.2, 5400.2, 5500.1, 5600.0,
                    5700.0, 5800.1, 5900.1, 6000.1, 6100.1, 6200.1, 6300.0, 6400.0, 6500.0, 6599.8,
                    6700.0, 6799.4, 6900.2, 7000.5, 7099.3, 7199.8, 7300.4, 7399.9, 7499.6, 7594.6,
                    7702.8, 7799.6, 7899.9, 7999.8, 8098.8, 8199.4, 8301.1, 8399.9, 8499.8, 8599.6,
                    8699.5, 8799.6, 8898.6, 8999.3, 9099.5, 9695.5])

r  = table[:, 2]
dr = table[:, 2 + 5]

nbin = 3
bins = np.linspace(r.min(), r.max(), nbin + 1)
bcen = (bins[1:] + bins[:- 1]) * 0.5

bin_ind = np.digitize(r, bins[:- 1])

plt.plot(SDSS[:, 0], SDSS[:, 1] / SDSS[:, 1].max() * 4, "--k", lw = 1)
plt.plot(JPAS[:, 0], JPAS[:, 1] / JPAS[:, 1].max() * 4, "-", color = "#FF4500", lw = 0.5)

output   = open("../../inputs/SDSS_sigma.txt", "w")
jpas_out = open("../../inputs/JPAS_sigma.txt", "w")
output.write("{0:12d}{1:12d}\n".format(nbin, wlength.size))
jpas_out.write("{0:12d}{1:12d}\n".format(nbin, wl_jpas.size))
output.write((wlength.size * "{:12.1f}" + "\n").format(*wlength))
jpas_out.write((wl_jpas.size * "{:12.1f}" + "\n").format(*wl_jpas))
for i in range(1, nbin + 1) :
	medians = np.median(sigmas[bin_ind == i], axis = 0)
	mean    = np.median(fluxes[bin_ind == i, 2] / sigmas[bin_ind == i, 2])

	output.write(("{:12.1f}" + medians.size * "{:12.4e}" + "\n").format(mean, *medians))

	plt.plot(wlength, medians / medians[2], "o-", label = str(round(mean/100.) * 100), mew = 0)

	polint = inter(wlength, medians, k = 2)
	tr_jpas = polint(wl_jpas)
	plt.plot(wl_jpas, tr_jpas / tr_jpas[np.argmin(abs(wl_jpas - wlength[2]))], "-k", lw = 2)

	jpas_out.write(("{:12.1f}" + tr_jpas.size * "{:12.4e}" + "\n").format(mean, *tr_jpas))

output.close()
jpas_out.close()

plt.gca().set_xticks(wlength)
plt.gca().set_xticklabels(["u'", "g'", "r'", "i'", "z'"])
plt.xlabel("wavelength [A]")
plt.ylabel("normalized flux at r'")
plt.legend(loc = 0, frameon = False, numpoints = 1, title = "~(S/N)_r'")
plt.tight_layout()
plt.show()
