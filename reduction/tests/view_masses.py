#!/usr/bin/python

import matplotlib.pyplot as plt
import numpy as np
import os

#sp = np.loadtxt("../inputs/real_deal3/sample_pars_m.log", dtype = np.int, usecols = (0, 1, 2))
#ja = np.loadtxt("../inputs/real_deal3/sample_pars_m.log", usecols = (25, 26, 27))
#rs = np.loadtxt("../inputs/real_deal3/sample_pars_m.log", usecols = (5,))
#mg = np.loadtxt("../inputs/real_deal3/sample_pars_m.log", usecols = (7, 9))
#ex = np.loadtxt("../inputs/real_deal3/sample_pars_m.log", usecols = (17, 19))
#pa = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_kcorr_blanton_model_ugriz")])
#kc = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_kcorr_model_ugriz")])
#
#sp = np.loadtxt("../inputs/real_deal3/sample_pars.log", dtype = np.int, usecols = (0, 1, 2))
#ja = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (25, 26, 27))
#rs = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (5,))
#mg = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (7, 9))
#pa = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_ugriz")])
#kc = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_kcorr_ugriz")])
#
sp = np.loadtxt("../inputs/real_deal3/sample_pars.log", dtype = np.int, usecols = (0, 1, 2))
ja = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (25, 26, 27))
rs = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (5,))
mg = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (7, 9))
ex = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (17, 19))
pa = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_kcorr_blanton_ugriz")])
kc = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.endswith(".log") and file.startswith("dynbasfit_ugriz")])

kcor = -2.5 * np.log10(np.loadtxt("../inputs/real_deal3/maggies.dat")[:, 1:] / np.loadtxt("../inputs/real_deal3/maggies_z0p00_at_z0p10_filters.dat")[:, 1:])[:, [0, 2]]
mg = mg - ex - kcor

def fu(x) :
  try :
    return int(x)
  except SyntaxError :
    return int(x[1:])

dyn_rmasses = []
dyn_kmasses = []
sel_names   = []
chisr       = []
chisk       = []

for file1, file2 in zip(pa, kc) :
  f1 = open(file1, "r")
  f2 = open(file2, "r")

  trip = map(fu, f1.readline()[:-1].split("=")[1].split("spSpec-")[1].split(".txt")[0].split("-"))
  trip = map(fu, f2.readline()[:-1].split("=")[1].split("spSpec-")[1].split(".txt")[0].split("-"))
  sel_names.append(trip)

  for i in xrange(15) : line1 = f1.readline()[:-1]
  for i in xrange(15) : line2 = f2.readline()[:-1]

  dyn_rmasses.append(map(eval, line1.split("=")[1].split()))
  dyn_kmasses.append(map(eval, line2.split("=")[1].split()))

  for i in xrange(5) : line1 = f1.readline()[:-1]
  for i in xrange(5) : line2 = f2.readline()[:-1]

  chisr.append(map(eval, line1.split("=")[1].split()))
  chisk.append(map(eval, line2.split("=")[1].split()))

chisr       = np.array(chisr)
chisk       = np.array(chisk)
sel_names   = np.array(sel_names)
mask        = [j for i in xrange(len(sel_names)) for j in xrange(sp.shape[0]) if np.all(sp[j] == sel_names[i])]
dyn_rmasses = np.array(dyn_rmasses)[np.array([chi == chi.min() for chi in chisr])]
dyn_kmasses = np.array(dyn_kmasses)[np.array([chi == chi.min() for chi in chisk])]
#dyn_rmasses = np.array(dyn_rmasses)[:, 0]
#dyn_kmasses = np.array(dyn_kmasses)[:, 0]
sp          = sp[mask]
jar_masses  = ja[mask]
rs          = rs[mask]
delta       = np.log10(dyn_rmasses) - jar_masses[:, 0]
ur          = (mg[:, 0] - mg[:, 1])[mask]
mask        = (jar_masses[:, 0] > 0.0)

plt.gca().set_aspect("equal")
plt.plot([8, 12], [8, 12], "--r", lw = 1)
plt.errorbar(np.log10(dyn_rmasses[mask]), jar_masses[mask, 0], np.array([jar_masses[mask, 0] - jar_masses[mask, 1], jar_masses[mask, 2] - jar_masses[mask, 0]]), fmt = "ko", ecolor = "k")
#plt.scatter(np.log10(dyn_rmasses[mask]), jar_masses[mask, 0], c = "k", lw = 0)
plt.scatter(np.log10(dyn_kmasses[mask]), jar_masses[mask, 0], c = "orange", lw = 0)

plt.xlabel(r"$\log M/M_\odot$ DynBaS $u'g'r'i'z'$", fontsize = 16)
plt.ylabel(r"$\log M/M_\odot$ by Brinchmann+2004", fontsize = 16)
plt.ylim(plt.xlim(8, 12))

#plt.figure()
#plt.plot(rs[mask], np.log10(dyn_rmasses[mask]) - jar_masses[mask, 0], "ok")
#plt.plot(rs[mask], np.log10(dyn_kmasses[mask]) - jar_masses[mask, 0], "o", color = "orange")
#plt.axvline(0.1, ls = "--", lw = 1.5, color = "r")
#plt.axhline(ls = "--", color = "r")
#plt.tight_layout()

plt.figure()
plt.hist([ur, ur[(abs(delta) <= 0.2)]], histtype = "barstacked")
#plt.hist(rs, histtype = "stepfilled", color = "k")
plt.tight_layout()

plt.show()
