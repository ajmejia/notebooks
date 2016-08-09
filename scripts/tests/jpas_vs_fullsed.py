import matplotlib.pyplot as plt
import numpy as np
import os

sp = np.loadtxt("../inputs/real_deal3/sample_pars.log", dtype = np.int, usecols = (0, 1, 2))
ja = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (25, 26, 27))
rs = np.loadtxt("../inputs/real_deal3/sample_pars.log", usecols = (5,))

jp = sorted([os.path.join(root, file) for root, subs, files in os.walk("jarle_consistency") for file in files if file.endswith(".log") and file.startswith("dynbasfit_processed_JPAS")])
fs = sorted([os.path.join(root, file) for root, subs, files in os.walk("jarle_consistency") for file in files if file.endswith(".log") and file.startswith("dynbasfit_processed_fullsed")])

def fu(x) :
  try :
    return int(x)
  except SyntaxError :
    return int(x[1:])

jp_masses, jp_mwage, jp_fwage, jp_met, jp_av = [], [], [], [], []
fs_masses, fs_mwage, fs_fwage, fs_met, fs_av = [], [], [], [], []
chis_jp = []
chis_fs = []

for fjp, ffs in zip(jp, fs) :
  fjp = open(fjp, "r")
  ffs = open(ffs, "r")

  for i in xrange(15) : linejp = fjp.readline()[:-1]
  for i in xrange(15) : linefs = ffs.readline()[:-1]

  jp_masses.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))
  jp_mwage.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))
  jp_fwage.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))
  jp_met.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))
  jp_av.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))

  fs_masses.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))
  fs_mwage.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))
  fs_fwage.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))
  fs_met.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))
  fs_av.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))

  chis_jp.append(map(eval, fjp.readline()[:-1].split("=")[1].split()))
  chis_fs.append(map(eval, ffs.readline()[:-1].split("=")[1].split()))

jp_masses = np.array(jp_masses)
jp_mwage = np.array(jp_mwage)
jp_fwage = np.array(jp_fwage)
jp_met = np.array(jp_met)
jp_av = np.array(jp_av)

fs_masses = np.array(fs_masses)
fs_mwage = np.array(fs_mwage)
fs_fwage = np.array(fs_fwage)
fs_met = np.array(fs_met)
fs_av = np.array(fs_av)

chis_jp   = np.array(chis_jp)
chis_fs   = np.array(chis_fs)

plt.figure()
plt.plot([6, 12], [6, 12], "--r", lw = 1)
plt.plot(np.log10(fs_masses[:, chis_fs == chis_fs.min(axis = 1)]), np.log10(jp_masses[:, chis_jp == chis_jp.min(axis = 1)]), "ok")
plt.xlabel(r"$\log M/M_\odot$ full SED")
plt.ylabel(r"$\log M/M_\odot$ J-PAS")
plt.tight_layout()

plt.figure()
plt.plot([6, 11], [6, 11], "--r", lw = 1)
plt.plot(fs_mwage[:, chis_fs == chis_fs.min(axis = 1)], jp_mwage[:, chis_jp == chis_jp.min(axis = 1)], "ok")
plt.xlabel(r"$<\log{t}>_M$ full SED")
plt.ylabel(r"$<\log{t}>_M$ J-PAS")
plt.tight_layout()

plt.figure()
plt.plot([6, 11], [6, 11], "--r", lw = 1)
plt.plot(fs_fwage[:, chis_fs == chis_fs.min(axis = 1)], jp_fwage[:, chis_jp == chis_jp.min(axis = 1)], "ok")
plt.xlabel(r"$<\log{t}>_F_r$ full SED")
plt.ylabel(r"$<\log{t}>_F_r$ J-PAS")
plt.tight_layout()

plt.figure()
plt.plot([0., 3], [0., 3], "--r", lw = 1)
plt.plot(fs_met[:, chis_fs == chis_fs.min(axis = 1)], jp_met[:, chis_jp == chis_jp.min(axis = 1)], "ok")
plt.xlabel(r"$<\log{Z/Z_\odot}>$ full SED")
plt.ylabel(r"$<\log{Z/Z_\odot}>$ J-PAS")
plt.tight_layout()

plt.show()
