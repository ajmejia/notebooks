
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
import pyfits as pyf
from string import strip

import os

set2_dir = "../inputs/SFHs_set2/"

sample = {}
for line in open(set2_dir + "ssag.log") :
	name_bc, name_ic = line[:-1].split("=")

	sample[strip(name_bc) + ".sed"] = "thread" + strip(name_ic) + ".fits.gz"

mass_g, mass_i = [], []
for j in xrange(10) :
	#fig, axs = plt.subplots(ncols = 2, nrows = 5, figsize = (7, 9))
	#axs = np.ravel(axs)
	for i, key in enumerate(sample.keys()[j*10:(j + 1)*10]) :

		s = open(set2_dir + "seds/" + key, "r")
		lines = s.readlines()
		
		f = pyf.open(set2_dir + sample[key])
		wl, flux = np.loadtxt(set2_dir + "seds/" + key, usecols = (0, 2), unpack = True)

		mass_g.append(eval(lines[33].split("=")[1][:-4]) * 1e9)
		mass_i.append(f[0].header["tform"] - f[0].header["burstage"])

		#axs[i].plot(f[1].data["wavelength"], f[1].data["flux"] / f[1].data["flux"][5000], "-r")
		#axs[i].plot(wl, flux / f[1].data["flux"][5000], "-b")
		#axs[i].text(4000, 2.5, key[:-4])
#
		#axs[i].set_xlim(3000, 9000)
		#axs[i].set_ylim(0, 3)
#
		#if i > 0  : axs[i].set_yticks(axs[i].get_yticks()[:-1])
		#if i == 8 : axs[i].set_xticks(axs[i].get_xticks()[:-1])
		#if i in [1, 3, 5, 7, 9] : axs[i].set_yticklabels([])
		#if i < 8 : axs[i].set_xticklabels([])
#
	#fig.text(0.45, 0.03, "wavelength")
	#axs[4].set_ylabel("normalized flux")
	#fig.subplots_adjust(wspace = 0, hspace = 0)
	#plt.savefig("sample_seds" + str(j) + ".png", bbox_inches = "tight")

mass_g, mass_i = np.array(mass_g), np.array(mass_i)

plt.figure()
plt.hist((mass_g - mass_i) / mass_i * 100., fc = "none", ec = "k", lw = 2)
plt.show()
