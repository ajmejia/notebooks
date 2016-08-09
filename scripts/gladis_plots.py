#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  gladis_plots.py
#  
#  Copyright 2014 Alfredo J. Mejia <mejia@cida.ve>
#  
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#  
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#  
#  

import numpy as np
import matplotlib.pyplot as plt

def main():

  def snr_dist():
    fig, axs = plt.subplots(1, 1, sharex=True, figsize=(8, 8))
  
    axs.hist(table["sn1_r"], 30, histtype="step", hatch="/", color="#7C45B9", range=(10, 60))
  
    axs.set_xlabel("SNR$r$", size=16)
    axs.set_ylabel("Frecuencia", size=16)

    plt.tight_layout()
    plt.savefig("SNR_dist.eps", bbox_inches="tight")
  
  def several_seds():
    from itertools import product
  
    fig, axs = plt.subplots(3, 2, sharex=True, sharey=True, figsize=(12, 8))
  
    for i, j in product(xrange(axs.shape[0]), xrange(axs.shape[1])):
      if j == 0: lb, =axs[i, j].plot(wls[i][j], fluxs[i][j], color="#5151CB", lw=1)
      else: lr, =axs[i, j].plot(wls[i][j], fluxs[i][j], color="#B04F4F", lw=1)
  
    axs[0, 0].text(0.95, 0.95, "SNR = 20", fontsize=14, ha="right", va="top", transform=axs[0, 0].transAxes)
    axs[1, 0].text(0.95, 0.95, "SNR = 40", fontsize=14, ha="right", va="top", transform=axs[1, 0].transAxes)
    axs[2, 0].text(0.95, 0.95, "SNR = 70", fontsize=14, ha="right", va="top", transform=axs[2, 0].transAxes)

    axs[0, 0].legend([lb, lr], ["Blue galaxy", "Red galaxy"], loc=2, prop=dict(size=14))
    axs[-1, 0].set_xlabel("Wavelength", size=16)
    axs[-1, 1].set_xlabel("Wavelength", size=16)
    axs[1, 0].set_ylabel("Normalized flux", size=16)

    plt.tight_layout()
    plt.savefig("several_seds.eps", bbox_inches="tight")
  
  def sigma_lambda():
    plt.figure()
    plt.plot(wl_s, sg, color="#4D4D4D", lw=1.5)
    plt.xlabel(r"Longitud de onda (\AA)", size=16)
    plt.ylabel("Sigma", size=16)
    plt.tight_layout()
    plt.savefig("sigma_lambda.eps", bbox_inches="tight")

  ddir = "../inputs/"

  table = np.genfromtxt(ddir + "snr.csv", dtype=None, names=True, delimiter=",")
  wl_s, sg = np.loadtxt("../inputs/SDSSERROR_median_filtered.DAT", usecols=(0,1), unpack=True)

  blues = ["processed_spSpec-51641-0314-227.txt",
           "processed_spSpec-51615-0303-201.txt",
           "processed_spSpec-51699-0349-587.txt"]
  reds = ["processed_spSpec-51671-0348-383.txt",
          "processed_spSpec-51660-0291-193.txt",
          "processed_spSpec-51662-0308-191.txt"]

  wls, fluxs = [], []
  for b, r in zip(blues, reds):
    wlb, flb = np.loadtxt("../inputs/SFGs/" + b, usecols=(0, 1), unpack=True)
    wlr, flr = np.loadtxt("../inputs/PGs/" + r, usecols=(0, 1), unpack=True)

    flb[(wlb>5420)&(wlb<5430)] = 0.
    flr[(wlr>5420)&(wlr<5430)] = 0.

    flb = flb / flb[np.argmin(np.abs(wlb - 6000))]
    flr = flr / flr[np.argmin(np.abs(wlr - 6000))]

    wls.append([wlb[flb>0], wlr[flr>0]])
    fluxs.append([flb[flb>0], flr[flr>0]])

  snr_dist()
  several_seds()
  sigma_lambda()

if __name__ == '__main__':
	main()

