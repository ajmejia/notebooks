import matplotlib.pyplot as plt
import pyfits as pyf
import numpy as np

names = np.loadtxt("../../inputs/SFHs_set3/set3_list.log", dtype = "|S")
d4000_sam, hdelta_sam = [], []
for name in names :
	f = pyf.open("../../inputs/SFHs_set3/" + name)

	d4000_sam.append(f[0].header["b4vn"])
	hdelta_sam.append(f[0].header["hdelta"])

d4000_lib, hdelta_lib = np.loadtxt("../../inputs/run09_sfh_catalog.txt", usecols = (45, 44), unpack = True)
d4000, hdelta, z = np.loadtxt("../../inputs//kauffmann/kauffmann.csv", delimiter = ",", skiprows = 1, usecols = (0, 1, 3), unpack = True)

zmask = z < 0.03

d4000_ssp = np.loadtxt("../../inputs/kauffmann/ssp.3color", usecols = (3,))
hdelta_ssp = np.loadtxt("../../inputs/kauffmann/ssp.7lsindx_sed_lick_system", usecols = (9,))
d4000_ssp1 = np.loadtxt("../../inputs/kauffmann/ssp1.3color", usecols = (3,))
hdelta_ssp1 = np.loadtxt("../../inputs/kauffmann/ssp1.7lsindx_sed_lick_system", usecols = (9,))
d4000_ssp2 = np.loadtxt("../../inputs/kauffmann/ssp2.3color", usecols = (3,))
hdelta_ssp2 = np.loadtxt("../../inputs/kauffmann/ssp2.7lsindx_sed_lick_system", usecols = (9,))
d4000_tau = np.loadtxt("../../inputs/kauffmann/tau.3color", usecols = (3,))
hdelta_tau = np.loadtxt("../../inputs/kauffmann/tau.7lsindx_sed_lick_system", usecols = (9,))
d4000_tau1 = np.loadtxt("../../inputs/kauffmann/tau1.3color", usecols = (3,))
hdelta_tau1 = np.loadtxt("../../inputs/kauffmann/tau1.7lsindx_sed_lick_system", usecols = (9,))
d4000_tau2 = np.loadtxt("../../inputs/kauffmann/tau2.3color", usecols = (3,))
hdelta_tau2 = np.loadtxt("../../inputs/kauffmann/tau2.7lsindx_sed_lick_system", usecols = (9,))

d4000_jar, hdelta_jar = np.loadtxt("../../inputs/kauffmann/jb_coadded.ori_colors", usecols = (4, 3), unpack = True)
d4000_mil, hdelta_mil = np.loadtxt("../../inputs/kauffmann/fort.15.43mil_gal", usecols = (4, 3), unpack = True)

plt.scatter(d4000[zmask], hdelta[zmask], s = 1, c = "#4D4D4D", alpha = 0.05)
plt.scatter(d4000_lib, hdelta_lib, s = 4, c = "#1E90FF", alpha = 0.1)
plt.plot(d4000_sam, hdelta_sam, ",", mec = "r", mew = 2)

plt.plot(d4000_mil, hdelta_mil, ",r")
plt.plot(d4000_jar, hdelta_jar, ".", color = "gray")

#plt.plot(d4000_ssp, hdelta_ssp, "-r", lw = 2)
#plt.plot(d4000_ssp1, hdelta_ssp1, "-g", lw = 2)
#plt.plot(d4000_ssp2, hdelta_ssp2, "-b", lw = 2)
#plt.plot(d4000_tau, hdelta_tau, "--r", lw = 2)
#plt.plot(d4000_tau1, hdelta_tau1, "--g", lw = 2)
#plt.plot(d4000_tau2, hdelta_tau2, "--b", lw = 2)

plt.xlim(0.7,2.5)
plt.ylim(-6,12)
plt.xlabel("Dn(4000)")
plt.ylabel("Hdelta")
plt.tight_layout()
plt.show()
