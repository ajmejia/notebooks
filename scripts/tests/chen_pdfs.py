#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st

labs = [r"$t_{\rm form}\,{\rm (Gyr)}$", r"$\gamma\,{\rm (1/Gyr)}$", r"$A$", r"$t_{\rm burst}\,{\rm (100\,Myr)}$", r"$\log(\tau)$", r"$Z/Z_\odot$", r"$\tau_V$", r"$\mu$"]

xs = [np.linspace(1.5, 13.5, 1000), np.linspace(0, 1, 1000), np.linspace(0.03, 4.0, 1000),
      np.linspace(0.3, 3.0, 1000), np.linspace(7, 9, 1000),
      np.linspace(0.02, 2.5, 1000), np.linspace(0, 6, 1000), np.linspace(0.1, 1, 1000)]
ys = [lambda x : np.ones(1000), lambda x : np.ones(1000), lambda x : - 0.66 ** x/(x * np.log(1 - 0.66)), lambda x : np.ones(1000),
      lambda x : np.ones(1000), lambda x : np.piecewise(x, [x < 0.2, x > 0.2], [0.05, 1.0]),
      lambda x : st.norm.pdf(x, 1.2, 0.98520478)*0.98520478*np.sqrt(2*np.pi)/2.19, lambda x : st.norm.pdf(x, 0.3, 0.36566883)*0.36566883*np.sqrt(2*np.pi)/0.62]

fig = plt.figure(figsize = (12, 4))

axs = [fig.add_subplot(2, 4, i + 1) for i in xrange(8)]

for i, ax in enumerate(axs) :
  ax.tick_params(labelleft = False, left = False, right = False)
  ax.set_xlim(xs[i].min(), xs[i].max())
  ax.set_ylim(0.0, 2.5)
  ax.set_xlabel(labs[i], fontsize = 16)

  if i == 2 : ax.set_xticks(np.linspace(0.0, 4.0, 5))
  if i == 7 : ax.set_xticks(np.linspace(0.2, 1.0, 5))
  if i == 8 : ax.set_xticks(np.linspace(50, 400, 6))

  func = ys[i]
  #ax.plot(xs[i], func(xs[i]), "-r", lw = 2)
  ax.fill_between(xs[i], func(xs[i])/func(xs[i]).max()+0.25, hatch = "//", facecolor = "none", edgecolor = "r")

fig.text(0.01, 0.5, r"${\rm PDF}$", fontsize = 14, rotation = "vertical")
fig.subplots_adjust(left = 0.04, bottom = 0.15, right = 0.98, top = 0.96, wspace = 0.12, hspace = 0.48)
plt.savefig("../plots/chen_pdfs.png", bbox_inches = "tight")
plt.show()
