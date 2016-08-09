#!/usr/bin/python

import pyfits as pyf
import numpy as np

path  = "./"
names = np.loadtxt(path + "fits_name_sample.log")

gr = []
for name in names :
	f = pyf.open(path + name)

	gr.append(f[0].header["mag_g"] - f[0].header["mag_r"])

np.savetxt("gr_fiber_color.txt", gr, fmt = "%12.5f", header = "(g - r)[mag]")
