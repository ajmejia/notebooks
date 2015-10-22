from matplotlib import rc

font = {"family":"sans-serif", "sans-serif":"Bitstream Vera Sans", "size":14.0}
text = {"usetex":False, "latex.preamble":r"\usepackage{amsmath},\usepackage{helvet},\renewcommand{\familydefault}{\sfdefault},\boldmath"}

rc("axes", linewidth=0.5, labelsize="medium", titlesize="medium", labelweight="bold", color_cycle="AFD65C D69D5C 5C8ED6 D6CE5C D67C5C 5CD6B5 D65CC1".split())
rc("xtick.major", width=0.5)
rc("xtick", labelsize="medium")
rc("ytick.major", width=0.5)
rc("ytick", labelsize="medium")
rc("legend", numpoints=1, fontsize="medium", frameon=False)
rc("lines", linewidth=1.0, markeredgewidth=0.0, markersize=7)
rc("patch", linewidth=0.0)
rc("font", **font)
rc("text", **text)
rc("savefig", dpi=92, bbox="tight")
