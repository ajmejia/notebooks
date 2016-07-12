from pylab import *
import pyfits as pyf

table_in = genfromtxt("set3+refill.log", dtype=None, names=True)
table_gs = genfromtxt("photoz3/true.parameters", dtype=None, names=True)

name_in = loadtxt("set3+refill.log", dtype=str, usecols=(0,))
name_gs = loadtxt("photoz3/true.parameters", dtype=str, usecols=(0,))
names = loadtxt("photoz3/set3.equiv", dtype=str, usecols=(0, 2))

sort_in = []
for name in names[:, 1]:
	for i, n in enumerate(name_in):
		if name == n:
			sort_in.append(i)
			break

sort_gs = []
for name in names[:, 0]:
	for i, n in enumerate(name_gs):
		if name == n:
			sort_gs.append(i)
			break

table_in = table_in[sort_in]
table_gs = table_gs[sort_gs]

name_in = name_in[sort_in]
name_gs = name_gs[sort_gs]

kws_in = "MASS MWLAGE RFWLA Z AV_EFF".split()
kws_gs = "Masszobs MWLAyr rFWLAyr ZZsun Av_eff".split()

#k=4
#
#if k in [1, 2]:
	#print name_in[(10**(table_gs[kws_gs[k]] - table_in[kws_in[k]]) - 1.0)*100. > 5.]
#else:
	#print name_in[abs(table_gs[kws_gs[k]] - table_in[kws_in[k]])/table_in[kws_in[k]]*100.0 > 5.]
#
#print table_in.shape, table_gs.shape
#plot(table_in[kws_in[k]], table_gs[[kws_gs[k]]], "ok")

fits = ['56643.595047', '56644.261118', '56646.331903', '56663.497440', '56642.459962']
equi = [47, 61, 94, 119, 129]
for i in xrange(5):
	if i < 4: f = pyf.open("SFHs_set3/" + fits[i] + ".fits.gz")
	else: pyf.open("SFHs_set4/thread2_" + fits[i] + ".fits.gz")
	
	#wl_sp, fl_sp = f[1].data["wavelength"], f[1].data["flux"]
	#wl_jp, fl_jp = loadtxt("photoz3/SSAG3" + "{0:03d}".format(equi[i]) + "_z0p00.sed", usecols=(0, 1), unpack=True)
	#
	#figure()
	#plot(wl_sp, fl_sp, "-k")
	#plot(wl_jp, fl_jp, "-r")
	#xlim(3000, 9000)
	#show()
	
	t_sp, sfr_sp = f[2].data["age"], f[2].data["sfr"]
	t_jp, sfr_jp = loadtxt("photoz3/SSAG3" + "{0:03d}".format(equi[i]) + "_z0p00.sfr", usecols=(0, 1), unpack=True)
	
	figure()
	plot(f[0].header["tform"] - t_sp, sfr_sp, "-k", label="SSAG")
	plot(f[0].header["tform"] - t_jp, sfr_jp, "-r", label="Gustavo")
	xlabel("look-back-time")
	ylabel("SFR")
	savefig("ejemplo.png", bbox_inches="tight")
	show()



