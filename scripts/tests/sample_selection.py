#!/usr/bin/python

from os import walk
from os.path import join as pjoin
import sys

from string import join, replace
import pyfits as pyf

#output_1 = open("fits_list.log", "w")
output_2 = open("test.dat", "w")

#keywords = ("RAOBJ",
            #"DECOBJ",
            #"SPEC1_G",
            #"SPEC1_R",
            #"SPEC1_I",
            #"SPEC2_G",
            #"SPEC2_R",
            #"SPEC2_I",
            #"VEL_DIS",
            #"VEL_DISE",
            #"Z",
            #"Z_ERR",
            #"Z_CONF")
keywords = ("RAOBJ", "DECOBJ")

fmt_str = join(["{0[" + str(i) + "]:>20s}" for i in xrange(len(keywords))])
output_2.write("#" + fmt_str.format(keywords) + "\n")

# Volume limited sample
zi, zf = 0.0, 1000.0

path_names = (pjoin(path, fits) for path, subdirs, names in walk(".") for fits in names if fits.endswith("fit.gz"))

# 1640320 fits files in DR7
ngalaxies   = 0
sample_size = 0

for nfits, fits_name in enumerate(path_names) :
  fits   = pyf.open(fits_name)
  if fits[0].header["OBJTYPE"] != "GALAXY" : continue
  ngalaxies += 1

  z = fits[0].header["Z"]

  if zi <= z <= zf :
    sample_size += 1

    cols = tuple([z if kw == "Z" else fits[0].header[kw] for kw in keywords])

    #output_1.write(fits_name + "\n")
    output_2.write(fmt_str.replace("s", ".14E").format(cols) + "\n")

#output_1.close()
output_2.close()

print "number of fits read      = {}".format(nfits + 1)
print "number of galaxies found = {}".format(ngalaxies)
print "sample size              = {}".format(sample_size)
