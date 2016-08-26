def AB2flux(ABmag, weff):
    return 10 ** (-0.4 * (ABmag - 22.5)) * 3.631e-6 * 1e-23 * (3e18 / weff ** 2)

def integrated_flux(SED, passband):
    mask = (passband[0, 0] <= SED[:, 0])&(SED[:, 0] <= passband[-1, 0])
    ipassband = interp(SED[mask, 0], passband[:, 0], passband[:, 1])

    return trapz(SED[mask, 1]*ipassband, SED[mask, 0])/trapz(ipassband, SED[mask, 0])
