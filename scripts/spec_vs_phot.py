#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as st
import thesis_redux_tools as rt

name_pho  = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", usecols = (0,), dtype = "|S")
name_sed  = np.loadtxt("../../outputs/spectroscopic_fit/names.log", dtype = "|S")
order_sed = [j for i in xrange(120) for j in xrange(12000) if name_pho[i] == name_sed[j][:12]+".fits.gz"]

tform, mass, burstm, burstm, hdelta, b4000 = np.loadtxt("../../inputs/set3_catalog.txt", usecols = (5, 16, 17, 18, 44, 45), unpack = True)

table = np.loadtxt("../../outputs/photometric_fit/remote_set3/photofit_SDSS.physical")
table[:, 2:6] = 10 ** table[:, 2:6]
table = rt.binned_stat(table, 100)

ave_sed = np.loadtxt("../../outputs/spectroscopic_fit/table_din.v3.log")[order_sed]
ave_sed[:, 1:3] = 10 ** ave_sed[:, 1:3]
ave_sed = rt.binned_stat(ave_sed, bin_size = 100)

#res_pho = [rt.err(table[:, i], table[:, i + 1]) if i in [0] else rt.err(table[:, i], table[:, i + 1], False) for i in xrange(0, 10, 2)]
#res_sed = [rt.err(table[:, ::2][:, i], ave_sed[:, i]) if i in [0] else rt.err(table[:, ::2][:, i], ave_sed[:, i], False) for i in xrange(5)]

# --------------------------------------------------------------------------------------------------
lm = np.array([0, 8])

plt.figure(figsize = (6.5, 6))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(ave_sed[:, 0], table[:, 1], s = 15, lw = 0, c = "k")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("spectroscopic mass")
plt.ylabel("photometric mass")
rt.percent_labels(plt.gca(), tks, tls)

plt.savefig("../../plots/photometric_fit/mass_phot_vs_spec.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------
#lm = np.array([8.5, 10.5])
lm = np.array([0, 14e9])

plt.figure(figsize = (6.5, 6))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(ave_sed[:, 1], table[:, 3], s = 15, lw = 0, c = "k")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("spectroscopic mass weighted log age")
plt.ylabel("photometric mass weighted log age")
rt.percent_labels(plt.gca(), tks, tls)

plt.savefig("../../plots/photometric_fit/mwage_phot_vs_spec.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------
#lm = np.array([8.5, 10.5])
lm = np.array([0, 14e9])

plt.figure(figsize = (6.5, 6))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(ave_sed[:, 2], table[:, 5], s = 15, lw = 0, c = "k")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("spectroscopic flux weighted log age")
plt.ylabel("photometric flux weighted log age")
rt.percent_labels(plt.gca(), tks, tls)

plt.savefig("../../plots/photometric_fit/fwage_phot_vs_spec.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------
lm = np.array([0.0, 2.5])

plt.figure(figsize = (6.5, 6))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(10 ** ave_sed[:, 3], table[:, 7], s = 15, lw = 0, c = "k")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("spectroscopic metallicity")
plt.ylabel("photometric metallicity")
rt.percent_labels(plt.gca(), tks, tls)

plt.savefig("../../plots/photometric_fit/metallicity_phot_vs_spec.png", bbox_inches = "tight")

# --------------------------------------------------------------------------------------------------
lm = np.array([0.0, 3.5])

plt.figure(figsize = (6.5, 6))
plt.plot(lm, lm, "--k")

tks, tls = [], []
for per in [10, 30, 50] :
  m1, m2 = rt.err_slope(per)

  tks.append(m2 * lm[1])
  tls.append((str(per) + r"\%"))

  plt.plot(lm, m1 * lm, "-.k")
  plt.plot(lm, m2 * lm, "-.k")

plt.scatter(ave_sed[:, 4], table[:, 9], s = 15, lw = 0, c = "k")

plt.xlim(lm)
plt.ylim(lm)
plt.xlabel("spectroscopic dust extinction")
plt.ylabel("photometric dust extinction")
rt.percent_labels(plt.gca(), tks, tls)

plt.savefig("../../plots/photometric_fit/extinction_phot_vs_spec.png", bbox_inches = "tight")

plt.show()
