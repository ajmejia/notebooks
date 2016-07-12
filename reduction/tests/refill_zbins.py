import numpy as np
import matplotlib.pyplot as plt

ptable = np.loadtxt("../../inputs/set3+red.log", usecols = (2,))
table4 = np.loadtxt("../../inputs/SFHs_set4/set4.log", usecols = (1,))
names  = np.loadtxt("../../inputs/SFHs_set4/set4.log", usecols = (0,), dtype = np.str)
lines  = np.array(open("../../inputs/SFHs_set4/set4.log", "r").readlines())[2:]

ind = range(table4.size)

mask = table4>12e9
ind  = np.random.choice(range(table4[mask].size), 20)
#print lines[0]

for line in lines[mask][ind] :
  print line[:-1]
