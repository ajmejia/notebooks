from matplotlib import rc

font = {"family":"serif", "serif":"Times New Roman", "size":9.0}
text = {"usetex":True, "latex.preamble":r"\usepackage{amsmath},\renewcommand{\familydefault}{\rmdefault},\renewcommand{\seriesdefault}{\bfseries},\boldmath"}

rc("axes", linewidth=0.5, labelsize=9.0, titlesize=9.0, labelweight="bold", color_cycle="AFD65C D69D5C 5C8ED6 D6CE5C D67C5C 5CD6B5 D65CC1".split())
rc("xtick.major", width=0.3)
rc("xtick", labelsize="small")
rc("ytick.major", width=0.3)
rc("ytick", labelsize="small")
rc("legend", numpoints=1, fontsize="xx-small", frameon=False)
rc("lines", linewidth=1.0, markeredgewidth=0.0, markersize=7)
rc("patch", linewidth=0.0)
rc("font", **font)
rc("text", **text)
rc("savefig", dpi=92, bbox="tight")

