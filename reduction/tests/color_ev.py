import matplotlib.pyplot as plt
import numpy as np

def f(m) :
	return 10 ** ( - 0.4 * (m + 21.10)) * 3111821.0187436366

t = 10 ** np.loadtxt("/home/alfredo/Dropbox/Modelos_Nidia/bc2003_Abr13/output.1ABmag", usecols = [0])
table = np.loadtxt("/home/alfredo/Dropbox/Modelos_Nidia/bc2003_Abr13/output.1ABmag", usecols = range(2, 7))

u = table[:, 0] + table[:, 1]
g = table[:, 0]
r = table[:, 0] - table[:, 2]
i = table[:, 0] - table[:, 3]
z = table[:, 0] - table[:, 4]

u = f(u)
g = f(g)
r = f(r)
i = f(i)
z = f(z)

du = np.diff(u)
dg = np.diff(g)
dr = np.diff(r)
di = np.diff(i)
dz = np.diff(z)

j = 90
k = 100

print t[j], t[k]

plt.plot([u[j], g[j], r[j], i[j], z[j]])
plt.plot([u[k], g[k], r[k], i[k], z[k]])
plt.show()
