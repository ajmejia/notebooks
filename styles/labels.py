params = dict(
    M=r"\bf $\log{M_*/\text{M}_\odot}$",
    log_t_M=r"\bf $\left<\log{t_*}/\text{yr}\right>_M$",
    log_t_L=r"\bf $\left<\log{t_*}/\text{yr}\right>_L$",
    log_Z_M=r"\bf $\left<\log{Z_*/\text{Z}_\odot}\right>_M$",
    log_Z_L=r"\bf $\left<\log{Z_*/\text{Z}_\odot}\right>_L$",
    Av=r"\bf $A_V$"
)
Delta = dict(
    M=r"\bf $\Delta\log{M_*}$",
    log_t_M=r"\bf $\Delta\left<\log{t_*}\right>_M$",
    log_t_L=r"\bf $\Delta\left<\log{t_*}\right>_L$",
    log_Z_M=r"\bf $\Delta\left<\log{Z_*}\right>_M$",
    log_Z_L=r"\bf $\Delta\left<\log{Z_*}\right>_L$",
    Av=r"\bf $\Delta A_V$"
)
Delta_int = dict((kw, Delta[kw].replace("Delta", r"Delta_\text{int}")) for kw in Delta)
delta = dict(
    M=r"\bf $\delta\log{M_*}$",
    log_t_M=r"\bf $\delta\left<\log{t_*}\right>_M$",
    log_t_L=r"\bf $\delta\left<\log{t_*}\right>_L$",
    log_Z_M=r"\bf $\delta\left<\log{Z_*}\right>_M$",
    log_Z_L=r"\bf $\delta\left<\log{Z_*}\right>_L$",
    Av=r"\bf $\delta A_V$"
)
delta_g05 = dict(
    M=r"\bf $\delta_{G05}\log{M_*}$",
    log_t_L=r"\bf $\delta_{G05}\left<\log{t_*}\right>_L$",
    log_Z_L=r"\bf $\delta_{G05}\left<\log{Z_*}\right>_L$"
)
