import matplotlib.pyplot as plt
import numpy as np

#data_1d = np.loadtxt("DINBAS1D.debug")
i, j, singular, successful = np.loadtxt("DINBAS2D.debug", unpack = True, dtype = "i4,i4,|S1,|S1")
#data_3d = np.loadtxt("DINBAS3D.debug")

singular   = singular == "T"
successful = successful == "T"

print "singular   {0}/{1}".format(singular.sum(), singular.size)
print "successful {0}/{1}".format(successful.sum(), successful.size)

plt.figure(figsize = (8, 8))
ax = plt.subplot(111, xlim = (1, 585), ylim = (1, 585), xlabel = "$i$", ylabel = "$j$")
ax.set_aspect("equal")
ax.grid()
ax.set_xticks(range(0, 585, 35))
ax.set_yticks(range(0, 585, 35))

ax.plot(i[~singular], j[~singular], ",k")
plt.savefig("non_singular1.png", bbox_inches = "tight")

#ax.plot(i[successful], j[successful], ",k")
#plt.savefig("successful1.png", bbox_inches = "tight")

#plt.show()
