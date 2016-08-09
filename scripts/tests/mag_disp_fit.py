
from numpy import *
from matplotlib.pyplot import *
from scipy.optimize import *
from thesis_redux_tools import binner
import sys, getopt

def residuals(p, y, x) :
	b, c = p
	return y - b * exp(c * x)
def peval(x, p) :
	b, c = p
	return b * exp(c * x)

# ASK FOR DATA IN FILE
try :
	filename, options = getopt.getopt(sys.getarv[1:], "fmt:h")
except getopt.error, msg :
	print "ERROR: {0}".format(msg)
	print __doc__
else :
	for o, a in options :
		if o == "-fmt" : form = a
		else           : print __doc__


# LOAD DATA
table  = np.loadtxt(filename)
nbands = table.shape

guess        = [0.01, 0.02]
fitting_pars = []
for j in range(nbands) :
	# COMPUTE BINS
	mag, err = binner(table[:, j], table[:, j + nbands], 100, "median", ebar = True)
	# COMPUTE FIT
	pars  = leastsq(residuals, guess, args = (err, mag))
	fitting_pars.append(pars)

# WRITE OUTPUT FILE
fmt = "#" + len(bands) * "%15s"
header = fmt%bands

fmt = fmt.replace("s", ".4f")
fmt = fmt.replace("#", " ")

np.savetxt(filename[:index(".")] + "_fitting_pars.txt", fitting_pars, fmt = fmt, header = header)
