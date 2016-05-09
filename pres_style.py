from matplotlib import rc
from cycler import cycler

clist  = "#114477 #117755 #E8601C #771111 #771144 #4477AA #44AA88 #F1932D #AA4477 #774411 #777711 #AA4455".split()
ccycle = cycler("color", clist)

rc("figure", figsize=(6, 6))
rc("text", usetex=False)
rc("font", **{"family":"sans-serif", "sans-serif":"Bitstream Vera Sans", "size":16, "weight":700})
rc("axes", linewidth=1.0, labelsize="medium", titlesize="medium", labelweight=700, prop_cycle=ccycle)
rc("xtick", labelsize="xx-small")
rc("ytick", labelsize="xx-small")
rc("lines", linewidth=2.0, markeredgewidth=0.0, markersize=7)
rc("patch", linewidth=0.0)
rc("legend", numpoints=1, scatterpoints=1, fontsize="xx-small", handletextpad=0.4, handlelength=1, handleheight=1, frameon=False)
rc("savefig", dpi=92, format="png")

