#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import thesis_redux_tools as rt
import sys, os

ntables = sorted([os.path.join(root, file) for root, subs, files in os.walk("../../outputs/photometric_fit/photoz3/") for file in files if file.startswith("photofit_z") and file.endswith(".physical")])

z        = np.append(np.arange(0.0, 1.0, 0.05), np.arange(1.0, 3.1, 0.1))
ptables  = []
residual = []

for name in ntables :

  table = np.loadtxt(name)
  table[:, 6:8] = np.log10(table[:, 6:8])

  ptables.append(table)

  res = [rt.err(table[:, i], table[:, i + 1]) if i in [0] else rt.err(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2)]
  res.pop(2)

  residual.append(res)

md = np.array(map(np.median, residual, [1] * len(residual)))

lines = plt.plot(z, md, "-")
plt.axhline(0.0, ls = "--", color = "k")

plt.legend(lines, ["mass", "age", "metallicity", "extinction"], loc = 0, prop = {"size" : 15})

plt.xlim(-0.05, 3.05)

plt.xlabel("Redshift" , fontsize = 15)
plt.ylabel("Residuals", fontsize = 15)

plt.tight_layout()
plt.savefig("../../plots/photometric_fit/residual_versus_z.pdf", bbox_inches = "tight")
plt.show()
