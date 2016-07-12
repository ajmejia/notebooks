#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

table = np.loadtxt("/home/alfredo/Programs/Maestria/input/hansson_sample/output.csv", delimiter = ",", skiprows = 1)

r = table[:, 17]
e = table[:, 22]
z = table[:, - 2]

mask = (z <= 0.03)
plt.plot(r, e, ",", color = "0.4", label = "galaxy sample")
plt.plot(r[mask], e[mask], ",", color = "r", label = "z < 0.03 sample")
xl = plt.xlim(12, 22)
yl = plt.ylim( - 0.001, 0.08)
plt.xlabel("fiberMag_r[mag]")
plt.ylabel("fiberMagErr_r[mag]")
#plt.savefig("fiber_r_photo_curv.png", bbox_inches = "tight")
plt.savefig("model_r_photo_curv.png", bbox_inches = "tight")
