params = dict(
    M=r"$\log{M_\star/\text{M}_\odot}$",
    log_t_M=r"$\left<\log{t_\star}/\text{a\~no}\right>_M$",
    log_t_L=r"$\left<\log{t_\star}/\text{a\~no}\right>_L$",
    log_Z_M=r"$\left<\log{Z_\star/\text{Z}_\odot}\right>_M$",
    log_Z_L=r"$\left<\log{Z_\star/\text{Z}_\odot}\right>_L$",
    Av=r"$A_V$"
)
Delta = dict(
    M=r"$\Delta\log{M_\star}$",
    log_t_M=r"$\Delta\left<\log{t_\star}\right>_M$",
    log_t_L=r"$\Delta\left<\log{t_\star}\right>_L$",
    log_Z_M=r"$\Delta\left<\log{Z_\star}\right>_M$",
    log_Z_L=r"$\Delta\left<\log{Z_\star}\right>_L$",
    Av=r"$\Delta A_V$"
)
Delta_int = dict((kw, Delta[kw].replace("Delta", r"Delta_\text{int}")) for kw in Delta)
delta = dict(
    M=r"$\delta\log{M_\star}$",
    log_t_M=r"$\delta\left<\log{t_\star}\right>_M$",
    log_t_L=r"$\delta\left<\log{t_\star}\right>_L$",
    log_Z_M=r"$\delta\left<\log{Z_\star}\right>_M$",
    log_Z_L=r"$\delta\left<\log{Z_\star}\right>_L$",
    Av=r"$\delta A_V$"
)
delta_g05 = dict(
    M=r"$\delta_\text{G05}\log{M_\star}$",
    log_t_L=r"$\delta_\text{G05}\left<\log{t_\star}\right>_L$",
    log_Z_L=r"$\delta_\text{G05}\left<\log{Z_\star}\right>_L$"
)

ssag_lb = r"SSAG"
sdss_lb = r"RB"
jpas_lb = r"RM"
spec_lb = r"RA"

jpaso_lb = r"RM$^*$"
speco_lb = r"RA$^*$"

sfgs_lb = r"GFEs"
pags_lb = r"GPas"
algs_lb = r"Todas"

tsize = "x-small"
