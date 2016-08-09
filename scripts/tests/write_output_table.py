#/usr/bin/python

import numpy as np
import os, sys

def fu(x) :
  try :
    return int(x)
  except SyntaxError :
    return int(x[1:])

#try :
  #fid = sys.argv[1]
#except IndexError :
  #print("usage: write_output_table.py files-ID")
  #sys.exit(1)

fid = "kcorr_blanton_ugriz"

file_list = sorted([os.path.join(root, file) for root, subs, files in os.walk(".") for file in files if file.startswith("dynbasfit_" + fid)])


sel_names, masses, m_ages, f_ages, metallicities, extinctions, chi_squares = [], [], [], [], [], [], []
for file_name in file_list :
  f = open(file_name, "r")

  trip = map(fu, f.readline()[:-1].split("=")[1].split("spSpec-")[1].split(".txt")[0].split("-"))
  sel_names.append(trip)

  for i in xrange(14) : line = f.readline()[:-1]

  masses.append(map(eval, f.readline()[:-1].split("=")[1].split()))
  m_ages.append(map(eval, f.readline()[:-1].split("=")[1].split()))
  f_ages.append(map(eval, f.readline()[:-1].split("=")[1].split()))
  metallicities.append(map(eval, f.readline()[:-1].split("=")[1].split()))
  extinctions.append(map(eval, f.readline()[:-1].split("=")[1].split()))
  chi_squares.append(map(eval, f.readline()[:-1].split("=")[1].split()))

  f.close()

header = "mjd plate fiber_ID mass_1d mass_2d mass_3d m_age_1d m_age_2d m_age_3d f_age_1d f_age_2d f_age_3d Z_1d Z_2d Z_3d Av_1d Av_2d Av_3d chi_1d chi_2d chi_3d".split()
header = ("#%5s %5s %8s" + 3 * "%12s" + 6 * "%10s" + 6 * "%7s" + 3 * "%15s" + "\n")%(tuple(header))

f = open(fid + "_dynbas_table.log", "w")
f.write(header)

table = np.hstack((sel_names, masses, m_ages, f_ages, metallicities, extinctions, chi_squares))

np.savetxt(f, table, fmt = " %05d  %04d      %03d" + 3 * "%12.3e" + 6 * "%10.6f" + 6 * "%7.2f" + 3 * "%15.4f")
f.close()

np.savetxt(fid + "_best_mass.log", table[:, [3, 4, 5]][np.array([chi == chi.min() for chi in np.array(chi_squares)])])
