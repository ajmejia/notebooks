#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np

c_label = ["u-g", "g-r", "r-i", "i-z"]
for i in xrange(4):
	F1 = str(120 + i)
	F2 = str(121 + i)
	fname = "cb2010_hr_milesx_o62_chab_ssp.color_F"+F1+"_F"+F2
	Z, col = np.loadtxt(fname, usecols = (0, 8), unpack = True)
	col_rf = col[0]
	residual = (col - col_rf) / col_rf * 100.0

	plt.plot(Z, residual, "-", label = c_label[i])

plt.legend(loc = 0, frameon = False)
xl = plt.xlim(0.0, 1)
plt.ylim(-100, 100)
plt.plot(xl, (0.0, 0.0), "--", color = "0.5")
plt.axvspan(0.02, 0.03, fc = "k", alpha = 0.1)
plt.xlabel(r"redshift", fontsize = 12)
plt.ylabel(r"residual (%)", fontsize = 12)
plt.show()
