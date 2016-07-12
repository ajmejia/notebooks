#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
from itertools import product

def max_diff(sed1, sed2) :
  #return np.abs((sed1 - sed2)).max()
  norm1 = np.sqrt(sum(sed1 * sed1))
  norm2 = np.sqrt(sum(sed2 * sed2))
  return np.arccos(np.sum(sed1 * sed2) / (norm1 * norm2)) * 180 / np.pi

table1 = np.loadtxt("../../codes/fort.90")
table2 = np.loadtxt("../../codes/fort.91")

t      = table1[:, 0]
#sed    = table1[:, 1:]
sed_op = table1[:, 876:6176]
sed_uv = table1[:, 1:6176]
sed_ir = table1[:, 876:]
pho    = table2[:, 1:]

x, y, diff_suv, diff_sop, diff_sir, diff_p = [], [], [], [], [], []
for i, j in product(xrange(t.size), xrange(t.size)) :
  x.append(i)
  y.append(j)
  diff_p.append(max_diff(pho[i], pho[j]))
  diff_sop.append(max_diff(sed_op[i], sed_op[j]))
  diff_suv.append(max_diff(sed_uv[i], sed_uv[j]))
  diff_sir.append(max_diff(sed_ir[i], sed_ir[j]))

x        = np.array(x)
y        = np.array(y)
diff_sop = np.array(diff_sop)
diff_suv = np.array(diff_suv)
diff_sir = np.array(diff_sir)
diff_p   = np.array(diff_p)

plt.figure()
plt.title("$u'g'r'i'z'$ photometry")
plt.scatter(x, y, c = diff_p, s = 7, lw = 0, vmin = min(diff_p), vmax = max(diff_p))
cb = plt.colorbar()

plt.xlabel("$i$ ages", fontsize = 14)
plt.ylabel("$j$ ages", fontsize = 14)
cb.set_label("$\max{|L_i - L_j|}$", fontsize = 14)
plt.xlim(0, 194)
plt.ylim(0, 194)
plt.savefig("SSPs_pho_comp.png", bbox_inches = "tight")

plt.figure()
plt.title("SDSS spectroscopy")
plt.scatter(x, y, c = diff_sop, s = 7, lw = 0, vmin = min(diff_p), vmax = max(diff_p))
cb = plt.colorbar()

plt.xlabel("$i$ ages", fontsize = 14)
plt.ylabel("$j$ ages", fontsize = 14)
cb.set_label("$\max{|L_i - L_j|}$", fontsize = 14)
plt.xlim(0, 194)
plt.ylim(0, 194)
plt.savefig("SSPs_spe_sdss_comp.png", bbox_inches = "tight")

plt.figure()
plt.title("UV+SDSS spectroscopy")
plt.scatter(x, y, c = diff_suv, s = 7, lw = 0, vmin = min(diff_p), vmax = max(diff_p))
cb = plt.colorbar()

plt.xlabel("$i$ ages", fontsize = 14)
plt.ylabel("$j$ ages", fontsize = 14)
cb.set_label("$\max{|L_i - L_j|}$", fontsize = 14)
plt.xlim(0, 194)
plt.ylim(0, 194)
plt.savefig("SSPs_spe_uvsdss_comp.png", bbox_inches = "tight")

plt.figure()
plt.title("SDSS+IR spectroscopy")
plt.scatter(x, y, c = diff_sir, s = 7, lw = 0, vmin = min(diff_p), vmax = max(diff_p))
cb = plt.colorbar()

plt.xlabel("$i$ ages", fontsize = 14)
plt.ylabel("$j$ ages", fontsize = 14)
cb.set_label("$\max{|L_i - L_j|}$", fontsize = 14)
plt.xlim(0, 194)
plt.ylim(0, 194)
plt.savefig("SSPs_spe_sdssir_comp.png", bbox_inches = "tight")

plt.show()
