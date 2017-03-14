from pylab import *
import itertools as it

def binner(x, y, bins=None, range=None, stats=(mean,lambda x: std(x,ddof=1))):
    if bins==None: bins = 10
    if range==None:
        mn, mx = x.min(), x.max()
    else:
        mn, mx = range

    ind = lexsort((y, x))
    xs, ys = x[ind], y[ind]

    xedges = linspace(mn, mx, bins+1)
    i_bins = digitize(xs, xedges[:-1])

    x_cen = (xedges[:-1]+xedges[1:]) * 0.5
    y_sta = array([stats[0](ys[i_bins==i]) for i in xrange(1, xedges.size)])

    myer = array([stats[1](ys[i_bins==i]) for i in xrange(1, xedges.size)])
    if len(myer.shape)==2:
        if myer.shape[0]==2:
            y_dev = myer
        elif myer.shape[1]==2:
            y_dev = myer.T
        else:
            raise ValueError, "y_dev has a shape {} incompatible with errorbar.".format(myer.shape)
    elif len(myer.shape)==1:
        y_dev = array([myer, myer])

    return x_cen, y_sta, y_dev

def last(nbin, span):
    return span[0], span[-1]+diff(span)/nbin

def nbins(sample, range=None) :
    IQR = lambda x: percentile(x, 75.0) - percentile(x, 25.0)
    if range==None:
        mn, mx = min(sample), max(sample)
    else:
        mn, mx = range

    mask = (sample>=mn)&(sample<=mx)
    binsize = (2 * IQR(sample[mask])/mask.sum() ** (1. / 3))

    return (mx-mn)/binsize, mn, mx, binsize

def err(xr, x, relative=True):
    '''
    Compute relative error between theoretical quantity x and real xr.
    '''

    if relative:
        return (x - xr)/xr
    if not relative:
        return x - xr

def redux_histogram2d(xbins, ybins, X, Y, Z, stat=average, fill=None):
    reduced_Z = zeros([xbins.size-1,ybins.size-1])
    for i, j in it.product(range(xbins.size-1), range(ybins.size-1)):
        mask = ((xbins[i]<X)&(X<=xbins[i+1]))&((ybins[j]<Y)&(Y<=ybins[j+1]))
        reduced_Z[i,j] = stat(Z[mask]) if mask.sum()>0 else np.nan
    if fill!=None: reduced_Z[isnan(reduced_Z)]=fill
    return reduced_Z
