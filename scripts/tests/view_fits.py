#!/usr/bin/python

import sys

if len(sys.argv[1:]) == 1 :
  file = sys.argv[1]
else :
  print "usage: view_fits dynbasfit_file"
  sys.exit(1)

try:
  with open(file, "r") as f :

    import numpy as np

    f.readline()
    f.readline()

    ages = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    mets = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    Av_s = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))

    f.readline()

    gen1d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    gen2d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    gen3d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    coe1d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    coe2d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    coe3d = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))

    f.readline()

    losvd = eval(f.readline()[:-1].split("=")[1])

    f.readline()

    masses  = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    mwlages = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    fwlages = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    mwlmets = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    Avs     = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))
    chis    = np.array(map(eval, f.readline()[:-1].split("=")[1].split()))

    wlength, f_obs, f_sig, f_m1d, f_m2d, f_m3d = np.loadtxt(f, unpack = True)
except IOError, msg :
  print msg
  sys.exit(1)

import matplotlib.pyplot as plt

fig = plt.figure()
ax  = fig.add_subplot(111)

ax.fill_between(wlength, f_obs - f_sig, f_obs + f_sig, linestyle = "dashed", lw = 1.5, facecolor = "gray", alpha = 0.7)
ax.plot(wlength, f_obs, "-x", color = "k",       lw = 1.5, label = "fitted SED")
ax.plot(wlength, f_m1d, "-x", color = "#DF6818", lw = 1.2, label = "DynBaS 1D ")
ax.plot(wlength, f_m2d, "-x", color = "#34A034", lw = 1.2, label = "DynBaS 2D ")
ax.plot(wlength, f_m3d, "-x", color = "#4B4B97", lw = 1.2, label = "DynBaS 3D ")
ax.legend(loc = 1, frameon = False)

ax.set_xlabel(r"$\lambda$\text{[\AA]}", fontsize = 16)
ax.set_ylabel(r"Flux\text{[L$\odot$/\AA]}", fontsize = 16)

plt.show()
