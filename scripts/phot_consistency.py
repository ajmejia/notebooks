#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as st
import thesis_redux_tools as rt

table = np.loadtxt("../../outputs/photometric_fit/photoz/photofit_JPASz3p00.physical")
table = rt.binned_stat(table, 50)

table[:, 2:6] = 10 ** table[:, 2:6]

res_pho = [rt.err(table[:, i], table[:, i + 1]) if i in [0] else rt.err(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2)]

#print table[np.argmax(res_pho[0]), 0], np.max(res_pho[0])
#print table[np.argmax(res_pho[1]), 2], np.max(res_pho[1])

mask = table[:, 6] > 0.4

# --------------------------------------------------------------------------------------------------
lm = np.array([0, 1])

plt.figure(figsize = (8, 7))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(table[mask, 0], table[mask, 1], s = 40, lw = 0, c = "b")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("input mass")
plt.ylabel("DynBaS mass")
rt.percent_labels(plt.gca(), tks, tls)

# --------------------------------------------------------------------------------------------------
lm = np.array([0, 1e9])

plt.figure(figsize = (8, 7))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(table[mask, 2], table[mask, 3], s = 40, lw = 0, c = "b")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("input mass weighted age")
plt.ylabel("DynBaS mass weighted age")
rt.percent_labels(plt.gca(), tks, tls)

# --------------------------------------------------------------------------------------------------
lm = np.array([0, 5e7])

plt.figure(figsize = (8, 7))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(table[mask, 4], table[mask, 5], s = 40, lw = 0, c = "b")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("input flux weighted age")
plt.ylabel("DynBaS flux weighted age")
rt.percent_labels(plt.gca(), tks, tls)

# --------------------------------------------------------------------------------------------------
lm = np.array([0.0, 2.5])

plt.figure(figsize = (8, 7))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(table[mask, 6], table[mask, 7], s = 40, lw = 0, c = "b")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("input metallicity")
plt.ylabel("DynBaS metallicity")
rt.percent_labels(plt.gca(), tks, tls)

# --------------------------------------------------------------------------------------------------
lm = np.array([0.0, 1.7])

plt.figure(figsize = (8, 7))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(table[mask, 8], table[mask, 9], s = 40, lw = 0, c = "b")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("input dust extinction")
plt.ylabel("DynBaS dust extinction")
rt.percent_labels(plt.gca(), tks, tls)

plt.show()
