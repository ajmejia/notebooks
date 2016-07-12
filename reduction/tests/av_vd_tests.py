#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np

ax = plt.subplot(111, xlim = (40, 400), ylim = (40, 400), xlabel = r"$\sigma_v$ simulado",
                                                          ylabel = r"$\sigma_v$ recuperado")
sub = ["CaII", "Mgb", "NaD", "SED"]
col = ["#1E90FF", "#FFA500", "#97657A", "#38942C"]

for k, s in enumerate(sub[:]) :

	SNs, tau, losvd_obs, losvd_ret = np.loadtxt("../../outputs/tests/losvd_test/dusty/losvd_test_" + s + ".txt", unpack = True)

	mask = (SNs == 10.0)

	x = np.array([np.mean(losvd_obs[mask][i * 10:(i + 1) * 10]) for i in xrange(23)])
	y = np.array([np.mean(losvd_ret[mask][i * 10:(i + 1) * 10]) for i in xrange(23)])
	i = np.lexsort((y, x))
	ax.plot(x[i], y[i], "o-", color = col[k], lw = 0.5, label = s)

	#x = np.array([np.mean(losvd_obs[mask][i * 10:(i + 1) * 10]) for i in xrange(23)])
	#y = np.array([np.mean(losvd_ret[mask][i * 10:(i + 1) * 10]) for i in xrange(23)])
	#i = np.lexsort((y, x))
	#ax.plot(x[i], y[i], "--", color = col[k], lw = 1.5)

ax.plot([40, 400], [40, 400], "--k", alpha = 0.5)
plt.legend(loc = 0, title = "spectral feature", numpoints = 1, frameon = False)
#plt.savefig("dusty_plain_vd_ret.pdf")
plt.show()
