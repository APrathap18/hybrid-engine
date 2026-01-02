import math
import numpy as np
from rocketcea.cea_obj_w_units import CEA_Obj

PSIA_TO_PA = 6894.757293168
POISE_TO_PAS = 0.1
MCAL_CM_S_K_TO_W_M_K = 0.4184  # RocketCEA default k unit -> W/(m·K)

# -----------------------------
# Compressible helpers
# -----------------------------
def area_mach(M, gamma):
    term = (2.0/(gamma+1.0)) * (1.0 + 0.5*(gamma-1.0)*M*M)
    exp  = (gamma+1.0) / (2.0*(gamma-1.0))
    return (1.0/M) * (term**exp)

def mach_from_area_ratio(A_ratio, gamma, supersonic=True, iters=80):
    if supersonic:
        lo, hi = 1.0001, 50.0
    else:
        lo, hi = 1e-6, 0.9999

    flo = area_mach(lo, gamma) - A_ratio
    fhi = area_mach(hi, gamma) - A_ratio
    # (Very rarely) if not bracketed, widen hi
    if supersonic and flo * fhi > 0:
        for _ in range(20):
            hi *= 2.0
            fhi = area_mach(hi, gamma) - A_ratio
            if flo * fhi <= 0:
                break
        else:
            raise RuntimeError("Could not bracket supersonic Mach solution for area ratio.")

    for _ in range(iters):
        mid = 0.5*(lo+hi)
        fmid = area_mach(mid, gamma) - A_ratio
        if fmid == 0:
            return mid
        if flo * fmid > 0:
            lo = mid
            flo = fmid
        else:
            hi = mid
    return 0.5*(lo+hi)

def p_over_p0(gamma, M):
    return (1.0 + 0.5*(gamma-1.0)*M*M)**(-gamma/(gamma-1.0))

def p_star_over_p0(gamma):
    return (2.0/(gamma+1.0))**(gamma/(gamma-1.0))

def Taw_from_T0(T0, gamma, Pr, M, turbulent=True):
    # recovery factor
    r = Pr**(1.0/3.0) if turbulent else Pr**0.5
    a = 0.5*(gamma-1.0)*M*M
    return T0 * (1.0 + r*a) / (1.0 + a)

# -----------------------------
# Updated wall_thickness
# -----------------------------
def wall_thickness(pc, of, eps, cstar, d_t, d_c, Tc, burn_time,
                   therm_cond, density, c_p, T_table, yield_strength,
                   oxname, fuelname, FS, chamb_length,
                   p_amb=101325.0, t_max0=0.05):
    """
    Returns:
      best: dict station->required thickness [m]
      h_g: dict station->h [W/m^2-K]
      stations: gas property dicts per station (as before)
      station_flow: dict station->{gamma, M, p_static, Taw}
    """

    cea = CEA_Obj(
        oxName=oxname,
        fuelName=fuelname,
        specific_heat_units="kJ/kg-K",
        viscosity_units="poise",
    )

    # --- Gas transport props from CEA ---
    Cp_c, mu_c, k_c, Pr_c = cea.get_Chamber_Transport(Pc=pc, MR=of)
    Cp_t, mu_t, k_t, Pr_t = cea.get_Throat_Transport(Pc=pc, MR=of)
    Cp_e, mu_e, k_e, Pr_e = cea.get_Exit_Transport(Pc=pc, MR=of, eps=eps)

    # Convert to SI
    Cp_c *= 1000.0; Cp_t *= 1000.0; Cp_e *= 1000.0
    mu_c *= POISE_TO_PAS; mu_t *= POISE_TO_PAS; mu_e *= POISE_TO_PAS
    k_c  *= MCAL_CM_S_K_TO_W_M_K
    k_t  *= MCAL_CM_S_K_TO_W_M_K
    k_e  *= MCAL_CM_S_K_TO_W_M_K

    d_e = d_t * math.sqrt(eps)

    stations = {
        "c": {"Cp": Cp_c, "mu": mu_c, "k": k_c, "Pr": Pr_c, "d": d_c},
        "t": {"Cp": Cp_t, "mu": mu_t, "k": k_t, "Pr": Pr_t, "d": d_t},
        "e": {"Cp": Cp_e, "mu": mu_e, "k": k_e, "Pr": Pr_e, "d": d_e},
    }

    # --- Get gamma from RocketCEA (fallback if method not available) ---
    # Prefer throat gamma as "constant gamma" for nozzle relations.
    try:
        _, gamma_c = cea.get_Chamber_MolWt_gamma(Pc=pc, MR=of)
    except Exception:
        gamma_c = 1.20

    try:
        _, gamma_t = cea.get_Throat_MolWt_gamma(Pc=pc, MR=of)
    except Exception:
        gamma_t = 1.20

    try:
        _, gamma_e = cea.get_Exit_MolWt_gamma(Pc=pc, MR=of, eps=eps)
    except Exception:
        gamma_e = gamma_t

    # --- Mach per station ---
    M_c = 0.0
    M_t = 1.0
    # Use gamma_t in area-Mach inversion (common practice for first-pass)
    M_e = mach_from_area_ratio(eps, gamma_t, supersonic=True)

    # --- Static pressure per station (for stress) ---
    pc_Pa = pc * PSIA_TO_PA
    p_c = pc_Pa
    p_t = pc_Pa * p_star_over_p0(gamma_t)
    p_e = pc_Pa * p_over_p0(gamma_t, M_e)

    # --- Taw per station (drive temp for convection) ---
    # Treat Tc as stagnation temperature T0 (first-pass assumption).
    T0 = float(Tc)
    Taw_c = Taw_from_T0(T0, gamma_c, Pr_c, M_c, turbulent=True)
    Taw_t = Taw_from_T0(T0, gamma_t, Pr_t, M_t, turbulent=True)
    Taw_e = Taw_from_T0(T0, gamma_e, Pr_e, M_e, turbulent=True)

    station_flow = {
        "c": {"gamma": gamma_c, "M": M_c, "p_static": p_c, "Taw": Taw_c},
        "t": {"gamma": gamma_t, "M": M_t, "p_static": p_t, "Taw": Taw_t},
        "e": {"gamma": gamma_e, "M": M_e, "p_static": p_e, "Taw": Taw_e},
    }

    # --- Gas-side h_g correlation ---
    # Use throat reference props for all stations (standard-ish Bartz-style usage)
    mu_ref, Cp_ref, Pr_ref = mu_t, Cp_t, Pr_t

    h_g = {}
    for s, p in stations.items():
        D = p["d"]
        h_g[s] = (
            (0.026 / (d_t ** 0.2))
            * ((mu_ref ** 0.2) * Cp_ref / (Pr_ref ** 0.6))
            * ((pc_Pa / cstar) ** 0.8)
            * ((d_t / D) ** 1.8)
        )

    # --- Thickness solve per station with bracketing expansion (Step 6) ---
    best = {}
    for s, p in stations.items():
        r_i_station = 0.5 * p["d"]

        # Use real h_g as the peak for dt suggestion
        dt = suggest_dt_explicit(
            r_i=r_i_station,
            t_wall=1e-4,
            L=chamb_length,
            rho_s=density,
            cp_s=c_p,
            k_s=therm_cond,
            n=40,
            hg_peak=h_g[s],
            safety=0.2,
        )

        # Drive convection with Taw and constant h
        Taw_s = station_flow[s]["Taw"]
        h_s   = h_g[s]
        Taw_func = (lambda t, Taw=Taw_s: Taw)
        h_func   = (lambda t, h=h_s: h)

        # Use station static pressure for stress (Step 3)
        p_i_station = station_flow[s]["p_static"]

        # Bracket-expand hi until survival
        hi = float(t_max0)
        for _ in range(20):
            t_vec, _, _, _, Twi, _ = simulate_wall_transient_cyl(
                r_i=r_i_station, t_wall=hi, L=chamb_length,
                rho_s=density, cp_s=c_p, k_s=therm_cond,
                t_end=burn_time, dt=dt,
                Tg_func=Taw_func, hg_func=h_func,
                outer_bc="adiabatic", To_func=None, ho_func=None,
                T_init=300.0, n=40
            )
            t_fail, _, _ = failure_time_for_thickness(
                p_i=p_i_station, p_o=p_amb, r_i=r_i_station, t_wall=hi,
                T_table=T_table, sy_table=yield_strength,
                FS=FS, t_vec=t_vec, Twi_hist=Twi
            )
            if t_fail is None:
                break
            hi *= 2.0
            if hi > 0.5:
                raise RuntimeError(f"[{s}] No survivable thickness found up to 0.5 m")

        # Now binary search within [t_min, hi]
        best[s] = thickness_for_target_time(
            p_i_station, p_amb, r_i_station,
            T_table, yield_strength,
            FS,
            chamb_length, density, c_p, therm_cond,
            burn_time, dt,
            Taw_func, h_func,
            outer_bc="adiabatic", To_func=None, ho_func=None,
            T_init=300.0, n=40,
            t_min=1e-4, t_max=hi, iters=30
        )

    # Step 7: return best + diagnostics
    return best, h_g, stations, station_flow

# ----------------------------
# Find dt
# ----------------------------
def suggest_dt_explicit(r_i, t_wall, L, rho_s, cp_s, k_s, n, hg_peak, safety=0.2):
    r_o = r_i + t_wall
    r_b = np.linspace(r_i, r_o, n + 1)
    r_c = 0.5 * (r_b[:-1] + r_b[1:])
    A_b = 2.0 * math.pi * r_b * L
    V   = math.pi * (r_b[1:]**2 - r_b[:-1]**2) * L

    # Diffusion limit
    alpha = k_s / (rho_s * cp_s)
    dr_min = np.min(np.diff(r_c))
    dt_diff = safety * (dr_min**2) / alpha

    # Inner boundary limit (cell 0)
    A_inner = A_b[0]
    A_int01 = A_b[1]
    dR01 = r_c[1] - r_c[0]
    G0 = k_s * A_int01 / dR01   # W/K
    dt_bc = safety * (rho_s * cp_s * V[0]) / (hg_peak * A_inner + G0)

    return min(dt_diff, dt_bc)


# # Find T_g and h_g functions
# def Tg_func(t, T0, T_hot, tau):
#     """
#     Docstring for Tg_func
    
#     :param t: Current time at which you want Tg
#     :param T0: Initial temp at t = 0
#     :param T_hot: Asymptotic final temperature
#     :param tau: Time constant (when you're 63% of the way from T_0 to T_hot)
#     :return: T_g
#     """
#     return T0 + (T_hot - T0) * (1.0 - math.exp(-t / tau))

# def hg_func(t, h0, h_hot, tau):
#     return h0 + (h_hot - h0) * (1.0 - math.exp(-t / tau))

# ----------------------------
# Yield strength interpolation
# ----------------------------
def yield_strength_at_T(T, T_table, sy_table):
    """
    T_table: temperatures [K]
    sy_table: yield strengths [Pa] at those temps
    Returns interpolated yield strength [Pa].
    """
    T_table = np.asarray(T_table, dtype=float)
    sy_table = np.asarray(sy_table, dtype=float)
    return float(np.interp(T, T_table, sy_table,
                           left=sy_table[0], right=sy_table[-1]))

# ----------------------------
# Thick-wall hoop stress (Lame)
# ----------------------------
def hoop_stress_inner(p_i, p_o, r_i, r_o):
    """
    Returns max hoop stress at inner radius for thick-walled cylinder
    using Lame's equation
    
    :param p_i: inner pressure (Pa)
    :param p_o: outer pressure (Pa)
    :param r_i: inner radius (m)
    :param r_o: outer radius (m)
    """

    if r_o <= r_i:
        raise ValueError("r_o must be greater than r_i")
    
    denom = (r_o**2 - r_i**2)

    A = (p_i * r_i**2 - p_o * r_o**2) / denom
    B = (r_i**2 * r_o**2 * (p_i - p_o)) / denom

    sigma_theta_i = A + B / (r_i**2)
    
    return sigma_theta_i

# --------------------------------------------
# Transient radial conduction in cylindrical wall
# --------------------------------------------
def simulate_wall_transient_cyl(r_i, t_wall, L, 
                                rho_s, cp_s, k_s, 
                                t_end, dt, 
                                Tg_func, hg_func, 
                                outer_bc = "adiabatic", 
                                To_func = None, ho_func = None, 
                                T_init = 300.0, n = 40
):
    """
    Simulates T(r,t) for a cylindrical wall from r_i to r_o=r_i+t_wall.

    Units (all SI):
      r_i, t_wall, L [m]
      rho_s [kg/m^3], cp_s [J/kg-K], k_s [W/m-K]
      t_end, dt [s]
      Tg_func(t)->[K], hg_func(t)->[W/m^2-K]
      outer convection: To_func(t)->[K], ho_func(t)->[W/m^2-K]
      Returns:
        t_vec [s], r_centers [m], T_hist [Nt x n], qin_hist [Nt], Twi_hist [Nt], Two_hist [Nt]
    """

    # Temperature varies only radially through the wall
    # The wall is a cylinder from inner radius r_i to outer radius r_o
    # Inner boundary heat input is convection from gas q_in(t) = h_g(T)(T_g(t)-T_w,i(t))
    # Outer boundary is adiabatic (no heat out) to solve for worst case scenario

    r_o = r_i + t_wall

    # n is the number of rings (annular control volumes) through the wall
    if n < 3:
        raise ValueError("Use n>=3 radial cells for stability/accuracy.")

    # Cell boundaries and centers
    # Separate the radii into small slices
    r_b = np.linspace(r_i, r_o, n + 1)               # boundaries
    r_c = 0.5 * (r_b[:-1] + r_b[1:])                 # centers
    dr_c = np.diff(r_c)                               # center-to-center spacing

    # Volumes and interface areas (cylindrical shell volumes)
    # Get volumes and areas at each cell radius
    V = math.pi * (r_b[1:]**2 - r_b[:-1]**2) * L      # cell volumes
    A_b = 2.0 * math.pi * r_b * L                     # boundary surface areas at each boundary radius

    Nt = int(math.floor(t_end / dt)) + 1              # Nt is number of time points saved (how many time steps to run)
    t_vec = np.linspace(0.0, dt * (Nt - 1), Nt)       # Create array of length Nt (timeline)

    T = np.full(n, float(T_init), dtype=float)        # Creates length-n array holding current temperature in each radial cell starting at T_init
    T_hist = np.zeros((Nt, n), dtype=float)           # Initialize array of n columns and Nt rows 

    # Initialize to zero
    qin_hist = np.zeros(Nt)   # inner heat flux [W/m^2] (positive into wall)
    Twi_hist = np.zeros(Nt)   # inner surface temperature [K] (cell 0)
    Two_hist = np.zeros(Nt)   # outer surface temperature [K] (cell n-1)

    # Precompute conduction "conductances" between adjacent cells:
    # q(i->i+1) = G[i] * (T[i] - T[i+1])
    # where G = k*A_interface / dR based off of Fourier's law for 1D heat conduction across a slab
    # Interface at r_b[i+1], distance between centers r_c[i+1]-r_c[i]
    # G is the thermal conductance in W/K between cell i and cell i + 1 (see equation for G above)
    # If T[i] > T[i+1], heat flows outward, otherwise heat flows inward
    # Precompute G
    G = np.zeros(n - 1, dtype=float)
    for i in range(n - 1):
        # Option A (your current approximation):
        A_int = A_b[i + 1]
        dR = r_c[i + 1] - r_c[i]
        G[i] = k_s * A_int / dR

        # Option B (more exact cylindrical form using centers):
        # G[i] = 2.0 * math.pi * k_s * L / math.log(r_c[i+1] / r_c[i])

    q_cond_plus = np.zeros(n, dtype=float)
    dT = np.zeros(n, dtype=float)

    for k in range(Nt):
        t = t_vec[k]
        T_hist[k, :] = T

        Tg = float(Tg_func(t))
        hg = float(hg_func(t))

        A_inner = A_b[0]
        q_conv_in = hg * A_inner * (Tg - T[0])
        qin_hist[k] = q_conv_in / A_inner
        Twi_hist[k] = T[0]

        A_outer = A_b[-1]
        if outer_bc == "adiabatic":
            q_conv_out = 0.0
        elif outer_bc == "convection":
            if To_func is None or ho_func is None:
                raise ValueError("To_func and ho_func required for outer convection BC.")
            To = float(To_func(t))
            ho = float(ho_func(t))
            q_conv_out = ho * A_outer * (To - T[-1])
        else:
            raise ValueError("outer_bc must be 'adiabatic' or 'convection'.")

        Two_hist[k] = T[-1]

        if k == Nt - 1:
            break  # don't advance past the final saved time

        # conduction
        q_cond_plus.fill(0.0)
        for i in range(n - 1):
            q_cond_plus[i] = G[i] * (T[i] - T[i + 1])

        # update
        dT.fill(0.0)

        net0 = q_conv_in - q_cond_plus[0]
        dT[0] = net0 * dt / (rho_s * cp_s * V[0])

        for i in range(1, n - 1):
            net = q_cond_plus[i - 1] - q_cond_plus[i]
            dT[i] = net * dt / (rho_s * cp_s * V[i])

        netN = q_cond_plus[n - 2] + q_conv_out
        dT[n - 1] = netN * dt / (rho_s * cp_s * V[n - 1])

        T = T + dT

    return t_vec, r_c, T_hist, qin_hist, Twi_hist, Two_hist

# ------------------------------------------------
# Burn-time / thickness check with temp-dependent yield
# ------------------------------------------------
def failure_time_for_thickness(
    p_i, p_o, r_i, t_wall,
    T_table, sy_table,
    FS,
    t_vec, Twi_hist
):
    """
    Checks when hoop stress at inner radius exceeds allowable (yield/FS),
    using inner wall temperature Twi_hist(t).
    Returns (t_fail, min_margin, margin_hist)
    """

    r_o = r_i + t_wall
    sigma_theta = hoop_stress_inner(p_i, p_o, r_i, r_o)          # Gets hoop stress at inner wall

    denom = (r_o**2 - r_i**2)
    sigma_z = (p_i * r_i**2 - p_o * r_o**2) / denom   # closed ends

    sigma_vm = math.sqrt(sigma_theta**2 + sigma_z**2 - sigma_theta*sigma_z)     # Axial stress (Von Mises)

    margin_hist = np.zeros_like(t_vec, dtype = float)           # Stores stress margin at every time point
    min_margin = float("inf")                                    # Sets min_margin to infinity so that the first computed margin will replace it and the worst-case margin will be tracked over time
    t_fail = None                                                # No failure found yet

    for i, t in enumerate(t_vec):
        sy = yield_strength_at_T(Twi_hist[i], T_table, sy_table) # [Pa], gets yield strength at the current inner wall temperature
        allow = sy / FS                                          # Allowable stress accounting for factor of safety
        margin = allow / sigma_vm                                # Available allowable stress / axial stress
        margin_hist[i] = margin                                  # If margin > 1 --> safe; margin = 1 --> right at allowable; margin < 1 --> failure
        min_margin = min(min_margin, margin)
        if t_fail is None and margin < 1.0:     
            t_fail = float(t)                                    # Record first failure time

    return t_fail, float(min_margin), margin_hist

# ------------------------------------------------
# Simple search: thickness needed to survive t_target
# ------------------------------------------------
def thickness_for_target_time(
    p_i, p_o, r_i,
    T_table, sy_table,
    FS,
    # thermal sim params:
    L, rho_s, cp_s, k_s,
    t_target, dt,
    Tg_func, hg_func,
    outer_bc="adiabatic", To_func=None, ho_func=None,
    T_init=300.0, n=40,
    t_min=1e-4, t_max=0.05, iters=30
):
    """
    Binary-search thickness [m] so that margin>=1 for all t in [0, t_target].
    """

    lo, hi = t_min, t_max

    # expand hi until it survives (or you hit a sanity cap)
    for _ in range(20):
        t_vec, _, _, _, Twi, _ = simulate_wall_transient_cyl(r_i=r_i, t_wall=hi, L = L, rho_s=rho_s, cp_s = cp_s, k_s = k_s, t_end = t_target, dt = dt, Tg_func = Tg_func, hg_func= hg_func, outer_bc="adiabatic", To_func = None, ho_func = None, T_init=300.0, n = 40)
        t_fail, _, _ = failure_time_for_thickness(p_i = p_i, p_o = p_o, r_i = r_i, t_wall=hi, t_vec=t_vec, Twi_hist=Twi)
        if t_fail is None:
            break
        hi *= 2.0
        if hi > 0.5:
            raise RuntimeError("No survivable thickness found up to 0.5 m")

    for _ in range(iters):
        mid = 0.5 * (lo + hi)

        t_vec, _, _, _, Twi, _ = simulate_wall_transient_cyl(r_i = r_i, t_wall = mid, L = L,
                                                             rho_s=rho_s, cp_s=cp_s, k_s=k_s,
                                                             t_end=t_target, dt = dt, 
                                                             Tg_func=Tg_func, hg_func = hg_func,
                                                             outer_bc=outer_bc, To_func=To_func, ho_func=ho_func,
                                                             T_init=T_init, n=n)
        
        t_fail, min_margin, _ = failure_time_for_thickness(
            p_i=p_i, p_o=p_o, r_i=r_i, t_wall=mid,
            T_table=T_table, sy_table=sy_table,
            FS=FS,
            t_vec=t_vec, Twi_hist=Twi
        )
        
        survives = (t_fail is None)
        if survives:
            best = mid
            hi = mid
        else:
            lo = mid

    
    if best is None:
        raise RuntimeError("Best thickness not set")

    return best