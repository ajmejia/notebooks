
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import pyfits as pyf

import os

set2_dir = "../inputs/SFHs_set2/"

weirds = {18:"thread11_56628.198931", 30:"thread9_56628.429931", 68:"thread14_56628.046314",
          81:"thread6_56627.830405", 83:"thread8_56629.201148", 86:"thread9_56628.283713",
          89:"thread10_56629.071848", 90:"thread5_56629.062899", 95:"thread9_56629.506443",
          96:"thread14_56628.367123"}

tau, mu, Av = [], [], []
fig, axs = plt.subplots(ncols = 2, nrows = 5, figsize = (7, 9))
axs = np.ravel(axs)
for i, key in enumerate(sorted(weirds)) :
	#s = open(set2_dir + "seds/SSAG00" + str(key) + ".sed", "r")
#
	#lines = s.readlines()
	#mass = eval(lines[3].split("=")[1])
	#tform = eval(lines[27].split("=")[1][:-4])

	f = pyf.open(set2_dir + weirds[key] + ".fits.gz")
	wl, flux = np.loadtxt(set2_dir + "seds/SSAG00" + str(key) + ".sed", usecols = (0, 2), unpack = True)

	print "%8.3f %8.3f %8.3f"%(f[0].header["A"]*f[0].header["mass"]/(1+f[0].header["A"]),f[0].header["burstm"], f[0].header["A"])

	tau.append(f[0].header["tau"])
	mu.append(f[0].header["mu"])
	Av.append(f[0].header["v"] - f[0].header["pv"])

	#print "%d %8.3f %8.3f"%(key, (mass - f[0].header["mass"]) / f[0].header["mass"] * 100, (tform - f[0].header["tform"]/1e9) / f[0].header["tform"]/1e9 * 100)

	#axs[i].plot(f[1].data["wavelength"], f[1].data["pflux"] / f[1].data["flux"][5000], "-k")
	axs[i].plot(f[1].data["wavelength"], f[1].data["flux"] / f[1].data["flux"][5000], "-r", label = r"$\tau = %.2f$ $\mu = %.2f$"%(tau[i], mu[i]))
	axs[i].plot(wl, flux / f[1].data["flux"][5000], "-b")
	#axs[i].legend(loc = 1, frameon = False)
	axs[i].text(4000, 3.5, str(key))

	axs[i].set_xlim(3000, 9000)
	axs[i].set_ylim(0, 4.5)

	if i > 0  : axs[i].set_yticks(axs[i].get_yticks()[:-1])
	if i == 8 : axs[i].set_xticks(axs[i].get_xticks()[:-1])
	if i in [1, 3, 5, 7, 9] : axs[i].set_yticklabels([])
	if i < 8 : axs[i].set_xticklabels([])

#print f[1].data["wavelength"][5000] = 7975

fig.text(0.45, 0.03, "wavelength")
axs[4].set_ylabel("normalized flux")

#fig, axs = plt.subplots(ncols = 1, nrows = 3)
#for i in xrange(len(sample)) :
	#f = pyf.open(sample[i])
	#
	#tau.append(f[0].header["tau"])
	#mu.append(f[0].header["mu"])
	#Av.append(f[0].header["v"] - f[0].header["pv"])

#axs[0].plot(tau[7:], mu[7:], "ok")
#axs[1].plot(tau[7:], Av[7:], "ok")
#axs[2].plot(mu[7:], Av[7:], "ok")
#
#axs[0].plot(tau[:7], mu[:7], "or")
#axs[1].plot(tau[:7], Av[:7], "or")
#axs[2].plot(mu[:7], Av[:7], "or")

#tau = np.array(tau)
#muu = np.array(mu)
#Av = np.array(Av)
#
#fig = plt.figure()

#plt.scatter(tau, mu, c = Av, lw = 0, s = 40)
#plt.xlabel(r"$\tau_V$", fontsize = 15)
#plt.ylabel(r"$\mu$", fontsize = 15)
#cb = plt.colorbar()
#cb.set_label(r"$A_V$", fontsize = 15)
#
#fig.tight_layout()
fig.subplots_adjust(wspace = 0, hspace = 0)
plt.show()
