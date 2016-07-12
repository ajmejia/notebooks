#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os, sys, getopt
import data_loader as dl

__doc__ = """Any help here would be hot!"""

try :
  opts_vals, args = getopt.getopt(sys.argv[1:], "h", longopts = ["redshift=", "help"])

  if args or not opts_vals :
    print __doc__
    sys.exit(1)

  for o, v in opts_vals :
    if   o == "-h" or o == "--help" :
      print __doc__
      sys.exit(0)
    elif o == "--redshift" :
      z = ("z{0:4.2f}".format(eval(v))).replace(".", "p")
except :
  print __doc__
  sys.exit(1)

f = open("../../inputs/photoz3/set3.counts")
z_id, counts = [], []
for line in f.readlines()[1:] :
  zi, ci = line[:-1].split()

  z_id.append(zi.replace(".", "p"))
  counts.append(eval(ci))

z_id   = np.array(z_id)
counts = np.array(counts)

jrun0 = dl.load_data(ID = "z0p00", n_SFH = counts[z_id == "z0p00"][0], n_trials = 50, n_bands = 56)
jrunz = dl.load_data(ID = z, n_SFH = counts[z_id == z][0], n_trials = 50, n_bands = 56)

itrialz = 103#np.random.randint(counts[z_id == z][0] * 50)
itrial0 = np.where(jrun0.library["SFH_name"] == jrunz.library["SFH_name"][itrialz].replace(z, "z0p00"))[0][0]

print
print "selected galaxy", itrialz, jrun0.library["SFH_name"][itrial0], jrunz.library["SFH_name"][itrialz]
print

def tra(val):
  if 220 <= val < 220*2:
    new = val - 220
  elif 220*2 <= val < 660:
    new = val - 220*2
  else:
    new = val

  return new

m3d0 = jrun0.dynbas[["coeff3d_1", "coeff3d_2", "coeff3d_3"]][itrial0]
it3d0 = map(tra, jrun0.dynbas[["gen3d_1", "gen3d_2", "gen3d_3"]][itrial0])
age3d0 = jrun0.dynbas[["age3d_1", "age3d_2", "age3d_3"]][itrial0]
wl_md0 = jrun0.models[["w_eff_"+str(i) for i in xrange(242,298)]][itrial0]
fl_md0 = jrun0.models[["L_mod_"+str(i) for i in xrange(242,298)]][itrial0]
wl_sy0 = jrun0.library[["w_eff_"+str(i) for i in xrange(242,298)]][itrial0]
fl_sy0 = jrun0.library[["L_syn_"+str(i) for i in xrange(242,298)]][itrial0]
tf0, tf_md0, tm0, tm_md0 = jrun0.physical["log_t_L_lib"][itrial0], jrun0.physical["log_t_L_mod"][itrial0], jrun0.physical["log_t_M_lib"][itrial0], jrun0.physical["log_t_M_mod"][itrial0]

m3dz = jrunz.dynbas[["coeff3d_1", "coeff3d_2", "coeff3d_3"]][itrialz]
it3dz = map(tra, jrunz.dynbas[["gen3d_1", "gen3d_2", "gen3d_3"]][itrialz])
age3dz = jrunz.dynbas[["age3d_1", "age3d_2", "age3d_3"]][itrialz]
wl_mdz = jrunz.models[["w_eff_"+str(i) for i in xrange(242,298)]][itrialz]
fl_mdz = jrunz.models[["L_mod_"+str(i) for i in xrange(242,298)]][itrialz]
wl_syz = jrunz.library[["w_eff_"+str(i) for i in xrange(242,298)]][itrialz]
fl_syz = jrunz.library[["L_syn_"+str(i) for i in xrange(242,298)]][itrialz]
tfz, tf_mdz, tmz, tm_mdz = jrunz.physical["log_t_L_lib"][itrialz], jrunz.physical["log_t_L_mod"][itrialz], jrunz.physical["log_t_M_lib"][itrialz], jrunz.physical["log_t_M_mod"][itrialz]

f = open("../../outputs/redshifted_ssps/fullsed_z0p00.log")

ages_fs  = np.array(map(eval, f.readline()[:-1].split()))
table_fs = np.loadtxt(f)

f.close()

f0 = open("../../outputs/redshifted_ssps/JPAS_z0p00.log")
fz = open("../../outputs/redshifted_ssps/JPAS_" + z + ".log")

ages_jp0  = np.array(map(eval, f0.readline()[:-1].split()))
table_jp0 = np.loadtxt(f0)
ages_jpz  = np.array(map(eval, fz.readline()[:-1].split()))
table_jpz = np.loadtxt(fz)

f0.close()
fz.close()

print "read {0} SSPs".format(table_fs.shape[1] - 1)

wl_fs   = table_fs[:, 0]
ssps_fs = table_fs[:, 1:]
wl_jp0   = table_jp0[:, 0]
ssps_jp0 = table_jp0[:, 1:]
wl_jpz   = table_jpz[:, 0]
ssps_jpz = table_jpz[:, 1:]

fig, axs = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12, 6))
axin0    = inset_axes(axs[0], width = "55%", height = "3%", loc = 1)
#axinz    = inset_axes(axs[1], width = "55%", height = "3%", loc = 1)

axs[0].set_xlim(3400, 9200)
axs[0].set_ylim(1e-8, 10.)

#print np.array(m3d0)/sum(np.array(m3d0))*100, np.array(m3dz), sum(m3dz)*100
#print it3d0, it3dz, ages_fs.size
#print np.log10(np.array(age3d0)), tm0, tm_md0, tf0, tf_md0, np.log10(np.array(age3dz)), tmz, tm_mdz, tfz, tf_mdz

mask0      = (ages_fs == age3d0[0])|(ages_fs == age3d0[1])|(ages_fs == age3d0[2])
maskz      = (ages_fs == age3dz[0])|(ages_fs == age3dz[1])|(ages_fs == age3dz[2])
lws       = np.ones(ssps_fs.shape[1] + 1) * 0.5
lws[-4:]  = 5

not_fitted_ssps0 = [list(zip(wl_jp0, ssps_jp0[:, i])) for i in xrange(ssps_fs.shape[1]) if ~mask0[i]]
fitted_ssps0     = [list(zip(wl_jp0, ssps_jp0[:, i])) for i in xrange(ssps_fs.shape[1]) if mask0[i]]
problem_sed0     = [list(zip(wl_sy0, fl_sy0))]

not_fitted_sspsz = [list(zip(wl_jpz, ssps_jpz[:, i])) for i in xrange(ssps_fs.shape[1]) if ~maskz[i]]
fitted_sspsz     = [list(zip(wl_jpz, ssps_jpz[:, i])) for i in xrange(ssps_fs.shape[1]) if maskz[i]]
problem_sedz     = [list(zip(wl_syz, fl_syz))]

lines0     = LineCollection(not_fitted_ssps0 + fitted_ssps0 + problem_sed0, linewidths = lws)
lines0.set_array(np.append(np.log10(ages_fs[~mask0]), np.append(np.log10(ages_fs[mask0]), tm0)))
linesz     = LineCollection(not_fitted_sspsz + fitted_sspsz + problem_sedz, linewidths = lws)
linesz.set_array(np.append(np.log10(ages_fs[~maskz]), np.append(np.log10(ages_fs[maskz]), tmz)))

wl_md0 = [item for item in wl_md0]
fl_md0 = [item for item in fl_md0]
wl_mdz = [item for item in wl_mdz]
fl_mdz = [item for item in fl_mdz]

axs[0].add_collection(lines0)
axs[0].plot(wl_md0, fl_md0, "-k", lw=0.8)
cb0 = plt.colorbar(lines0, cax = axin0, orientation="horizontal")
axs[1].add_collection(linesz)
axs[1].plot(wl_mdz, fl_mdz, "-k", lw=0.8)
#cbz = plt.colorbar(linesz, cax = axinz, orientation="horizontal")

axs[0].set_yscale("log")
axs[1].set_yscale("log")

axs[0].set_title("z0p00")
axs[0].set_xlabel(r"$\lambda[\text{\AA}]$",     fontsize = 16)
axs[0].set_ylabel(r"Flux[L$\odot/\text{\AA}]$", fontsize = 16)
cb0.set_label("Age[dex]", fontsize = 16)
axs[1].set_title("z3p00")
axs[1].set_xlabel(r"$\lambda[\text{\AA}]$",     fontsize = 16)
#cbz.set_label("Age[dex]", fontsize = 16)

plt.tight_layout()
plt.savefig("../../plots/draw_ssps_" + z + "_" + str(itrialz) + ".pdf", bbox_inches = "tight")
plt.show()
