#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import InterpolatedUnivariateSpline as inter

COSMOS = np.loadtxt("../../inputs/COSMOS_filters.txt")
JPAS   = np.loadtxt("../../inputs/J-PAS_filters.txt")
table  = np.genfromtxt("../../inputs/COSMOS/cosmos+scosmos.txt", missing = '""', filling_values = -9999.9)

mask_t = table[:, 3] == 0
mask_z = table[:, 2] < 0.03
mask_m = (table[:, 4] > - 99.)&(table[:, 5] > - 99.)&(table[:, 6] > - 99.)&(table[:, 7] > - 99.)&(table[:, 8] > - 99.)&(table[:, 9] > - 99.)&(table[:, 10] > - 99.)&(table[:, 11] > - 99.)&(table[:, 20] > - 99.)&(table[:, 22] > - 99.)
mask_V = (table[:, 6] > 12.)&(table[:, 6] < 25.)
mask   = (mask_m & mask_t & mask_V)
table  = table[mask]

wlength   = np.array([3881.6, 4345.3, 5401.8, 4619.4, 6165.2, 7591.2, 9153.1, 21337.8, 35075.1, 44365.8])
op_fluxes = 10 ** (- 0.4 * (table[:, 4:4+7] - 22.5)) * (3.631e-6 / wlength[:-3] ** 2) * (1e-23 * 3e18)
Ks_fluxes = 10 ** (- 0.4 * (table[:, 4+7] - 22.5)) * (3.631e-6 / wlength[7] ** 2) * (1e-23 * 3e18)
ir_fluxes = table[:, 20:20+4:2] * 1e-6 * 1e-23 * 3e18 / wlength[-2:] ** 2

op_sigmas = op_fluxes * table[:, 12:12+7] / (2.5 * np.log10(np.exp(1.0)))
Ks_sigmas = Ks_fluxes * table[:, 12+7] / (2.5 * np.log10(np.exp(1.0)))
ir_sigmas = table[:, [21,23]] * 1e-6 * 1e-23 * 3e18 / wlength[-2:] ** 2

fluxes = np.column_stack((op_fluxes, Ks_fluxes, ir_fluxes))
sigmas = np.column_stack((op_sigmas, Ks_sigmas, ir_sigmas))
wl_jpas = np.linspace(3485., 10075., 57)

V  = table[:, 6]
dV = table[:, 14]

nbin = 3
bins = np.linspace(V.min(), V.max(), nbin + 1)
bcen = (bins[1:] + bins[:- 1]) * 0.5

bin_ind = np.digitize(V, bins[:- 1])
fil_ind = np.argsort(wlength)

output = open("../../inputs/COSMOS_sigma.txt", "w")
output.write("{0:12d}{1:12d}\n".format(nbin, wlength.size))
output.write((wlength.size * "{:12.1f}" + "\n").format(*wlength))
for i in range(1, nbin + 1) :
	medians = np.median(sigmas[bin_ind == i], axis = 0)
	mean    = np.median(fluxes[bin_ind == i, 2] / sigmas[bin_ind == i, 2])

	output.write(("{:12.1f}" + medians.size * "{:12.4e}" + "\n").format(mean, *medians))

	plt.loglog(wlength[fil_ind], medians[fil_ind] / medians[fil_ind][2], "o-", label = str(round(mean/10.) * 10), mew = 0)
	#plt.plot(V[bin_ind == i], dV[bin_ind == i], ".")
	if i == 2 :
		polint = inter(wlength, medians, k = 2)
		tr_jpas = polint(wl_jpas)
		plt.plot(wl_jpas, tr_jpas / tr_jpas[17], "-k", lw = 2)

output.close()

plt.gca().set_xticks(wlength)
plt.gca().set_xticklabels(["u", "B", "g", "V", "r", "i", "z", "K", "3.6", "4.5"])
plt.xlabel("wavelength [A]")
plt.ylabel("sigma [erg/s/cm^2/A]")
plt.legend(loc = 0, frameon = False, numpoints = 1, title = "~(S/N)_V")
plt.tight_layout()
plt.show()
