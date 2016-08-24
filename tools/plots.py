from matplotlib import cm
import numpy as np
from numpy import ma

class MidpointNormalize(cm.colors.Normalize):
    # from http://matplotlib.org/devdocs/users/colormapnorms.html
    def __init__(self, vmin=None, vmax=None, midpoint=None, clip=False):
        self.midpoint = midpoint
        cm.colors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
        x, y = [self.vmin, self.midpoint, self.vmax], [0, 0.5, 1]
        return ma.masked_array(np.interp(value, x, y))
