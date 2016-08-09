#!/usr/bin/python

import pyfits as pyf
import os
import numpy as np

fits_list = (os.path.join(root, file) for root, subs, files in os.walk("run09/thread0/") for file in files if file.endswith(".fits.gz"))
tform = []
for fits in fits_list :
	f = pyf.open(fits)

	tform.append(f[0].header["tform"])

np.savetxt("SSAG_tform.log", tform, fmt = "%15.5e")
