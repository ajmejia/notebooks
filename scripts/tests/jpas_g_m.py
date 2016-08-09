
import matplotlib.pyplot as plt
import numpy as np
import os

gsols = sorted([os.path.join(root, file) for root, subs, files in os.walk("../../codes/.") for file in files if file.endswith(".nbp3")])
jsols = sorted([os.path.join(root, file) for root, subs, files in os.walk("../../codes/.") for file in files if file.startswith("dynbasfit_processed_JPAS")])
msols = sorted([os.path.join(root, file) for root, subs, files in os.walk("../../codes/.") for file in files if file.startswith("dynbasfit_processed_masked")])

ggen1, ggen2, ggen3, gcoe1, gcoe2, gcoe3, gchi1, gchi2, gchi3, gAv1, gAv2, gAv3 = [], [], [], [], [], [], [], [], [], [], [], []
n = []
for i in xrange(len(gsols)) :
  f = open(gsols[i], "r")
  for j in xrange(3) : f.readline()

  n.append(eval(f.readline()[:-1].split()[-1]))

  for j in xrange(8) : f.readline()

  db1 = np.array(f.readline()[:-1].split()[2:])

  gAv1.append(eval(db1[0]))
  gchi1.append(eval(db1[1]))
  gcoe1.append(eval(db1[-1]))
  ggen1.append(eval(db1[-2]))

  db2 = np.array(f.readline()[:-1].split()[2:])

  gAv2.append(eval(db2[0]))
  gchi2.append(eval(db2[1]))
  gcoe2.append(map(eval, db2[[-2, -1]]))
  ggen2.append(map(eval, db2[[-4, -3]]))

  db3 = np.array(f.readline()[:-1].split()[2:])

  gAv3.append(eval(db3[0]))
  gchi3.append(eval(db3[1]))
  gcoe3.append(map(eval, db3[[-3, -2, -1]]))
  ggen3.append(map(eval, db3[[-6, -5, -4]]))

gAv1  = np.array(gAv1 )
gchi1 = np.array(gchi1)
gcoe1 = np.array(gcoe1)
ggen1 = np.array(ggen1)

gAv2  = np.array(gAv2 )
gchi2 = np.array(gchi2)
gcoe2 = np.array(gcoe2)
ggen2 = np.array(ggen2)

gAv3  = np.array(gAv3 )
gchi3 = np.array(gchi3)
gcoe3 = np.array(gcoe3)
ggen3 = np.array(ggen3)

n     = np.array(n)

jgen1, jgen2, jgen3, jcoe1, jcoe2, jcoe3, jchi1, jchi2, jchi3, jAv1, jAv2, jAv3 = [], [], [], [], [], [], [], [], [], [], [], []
for i in xrange(len(gsols)) :
  f = open(jsols[i], "r")
  for j in xrange(6) : f.readline()

  jgen1.append(eval(f.readline()[:-1].split()[4:][0]))
  jgen2.append(map(eval, f.readline()[:-1].split()[4:]))
  jgen3.append(map(eval, f.readline()[:-1].split()[4:]))

  jcoe1.append(eval(f.readline()[:-1].split()[4:][0]))
  jcoe2.append(map(eval, f.readline()[:-1].split()[4:]))
  jcoe3.append(map(eval, f.readline()[:-1].split()[4:]))

  for j in xrange(7) : f.readline()
  Av = f.readline()[:-1].split()[3:]

  jAv1.append(eval(Av[0]))
  jAv2.append(eval(Av[1]))
  jAv3.append(eval(Av[2]))

  chi = f.readline()[:-1].split()[4:]

  jchi1.append(eval(chi[0]))
  jchi2.append(eval(chi[1]))
  jchi3.append(eval(chi[2]))

jAv1  = np.array(jAv1 )
jchi1 = np.array(jchi1)
jcoe1 = np.array(jcoe1)
jgen1 = np.array(jgen1)

jAv2  = np.array(jAv2 )
jchi2 = np.array(jchi2)
jcoe2 = np.array(jcoe2)
jgen2 = np.array(jgen2)

jAv3  = np.array(jAv3 )
jchi3 = np.array(jchi3)
jcoe3 = np.array(jcoe3)
jgen3 = np.array(jgen3)

mgen1, mgen2, mgen3, mcoe1, mcoe2, mcoe3, mchi1, mchi2, mchi3, mAv1, mAv2, mAv3 = [], [], [], [], [], [], [], [], [], [], [], []
for i in xrange(len(gsols)) :
  f = open(msols[i], "r")
  for j in xrange(6) : f.readline()

  mgen1.append(eval(f.readline()[:-1].split()[4:][0]))
  mgen2.append(map(eval, f.readline()[:-1].split()[4:]))
  mgen3.append(map(eval, f.readline()[:-1].split()[4:]))

  mcoe1.append(eval(f.readline()[:-1].split()[4:][0]))
  mcoe2.append(map(eval, f.readline()[:-1].split()[4:]))
  mcoe3.append(map(eval, f.readline()[:-1].split()[4:]))

  for j in xrange(7) : f.readline()
  Av = f.readline()[:-1].split()[3:]

  mAv1.append(eval(Av[0]))
  mAv2.append(eval(Av[1]))
  mAv3.append(eval(Av[2]))

  chi = f.readline()[:-1].split()[4:]

  mchi1.append(eval(chi[0]))
  mchi2.append(eval(chi[1]))
  mchi3.append(eval(chi[2]))

mAv1  = np.array(mAv1 )
mchi1 = np.array(mchi1)
mcoe1 = np.array(mcoe1)
mgen1 = np.array(mgen1)

mAv2  = np.array(mAv2 )
mchi2 = np.array(mchi2)
mcoe2 = np.array(mcoe2)
mgen2 = np.array(mgen2)

mAv3  = np.array(mAv3 )
mchi3 = np.array(mchi3)
mcoe3 = np.array(mcoe3)
mgen3 = np.array(mgen3)

#plt.plot(jcoe1, mcoe1, "rs", label = "1D")
#plt.plot(np.sum(jcoe2, axis = 1), np.sum(mcoe2, axis = 1), "g^", label = "2D")
#plt.plot(np.sum(jcoe3, axis = 1), np.sum(mcoe3, axis = 1), "ok", label = "3D")

#plt.plot(jgen1, mgen1, "rs", label = "1D")
#plt.plot(np.sum(jgen2, axis = 1), np.sum(mgen2, axis = 1), "g^", label = "2D")
#plt.plot(np.sum(jgen3, axis = 1), np.sum(mgen3, axis = 1), "ok", label = "3D")

plt.plot(gcoe1, jcoe1, "rs", label = "1D")
plt.plot(np.sum(gcoe2, axis = 1), np.sum(jcoe2, axis = 1), "g^", label = "2D")
plt.plot(np.sum(gcoe3, axis = 1), np.sum(jcoe3, axis = 1), "ok", label = "3D")
plt.plot([0, 3e11], [0, 3e11], "--r", lw = 1)
plt.xlabel(r"$M^\star$ de Gustavo", fontsize = 14)
plt.ylabel(r"$M^\star$ de Alfredo", fontsize = 14)

#plt.plot(gchi1, jchi1, "rs", label = "1D")
#plt.plot(gchi2, jchi2, "g^", label = "2D")
#plt.plot(gchi3, jchi3, "ok", label = "3D")
#plt.plot([0, 500], [0, 500], "--r", lw = 1)
#plt.xlabel(r"$\chi^2$ de Gustavo", fontsize = 14)
#plt.ylabel(r"$\chi^2$ de Alfredo", fontsize = 14)

plt.legend(loc = 0, frameon = False)
plt.tight_layout()
plt.show()
