import math
from rocketcea.cea_obj_w_units import CEA_Obj

PSIA_TO_PA = 6894.757293168
POISE_TO_PAS = 0.1
MCAL_CM_S_K_TO_W_M_K = 0.4184  # RocketCEA default k unit -> W/(m·K)

def wall_thickness(pc, of, eps, cstar, d_t, d_c, Tc, burn_time,
                   therm_cond, density, c_p, yield_strength,
                   oxname, fuelname):
    """
    pc = chamber pressure (psia)           [CEA expects psia]
    of = O/F ratio
    eps = expansion ratio Ae/At
    cstar = characteristic velocity (m/s) 
    d_t, d_c = throat/chamber diameter (m)
    Tc = combustion temperature (K)
    burn_time = burn time of engine (s)
    therm_cond = thermal conductivity of material (W/m-K)
    density = density of material (kg/m^3) 
    c_p = specific heat at constant pressure (J/kg-K) 
    yield_strength (Pa)
    """

    cea = CEA_Obj(
        oxName=oxname,
        fuelName=fuelname,
        specific_heat_units="kJ/kg-K",   # Cp returned in kJ/(kg·K)
        viscosity_units="poise",         # mu returned in poise
    )

    # --- Gas transport from CEA at stations ---
    # Returns: (Cp, mu (dynamic gas viscosity), k (thermal conductivity), Pr (prandtel's number))
    Cp_c, mu_c, k_c, Pr_c = cea.get_Chamber_Transport(Pc=pc, MR=of)
    Cp_t, mu_t, k_t, Pr_t = cea.get_Throat_Transport(Pc=pc, MR=of)
    Cp_e, mu_e, k_e, Pr_e = cea.get_Exit_Transport(Pc=pc, MR=of, eps=eps)

    # --- Convert to SI (Cp, mu). k conversion assumes RocketCEA default k units ---
    Cp_c = Cp_c * 1000.0          # kJ/kg-K -> J/kg-K
    Cp_t = Cp_t * 1000.0
    Cp_e = Cp_e * 1000.0

    mu_c = mu_c * POISE_TO_PAS    # poise -> Pa·s
    mu_t = mu_t * POISE_TO_PAS
    mu_e = mu_e * POISE_TO_PAS

    # If you need gas thermal conductivity in SI:
    k_c = k_c * MCAL_CM_S_K_TO_W_M_K  # -> W/(m·K) (only correct if k came in mcal/(cm·s·K))
    k_t = k_t * MCAL_CM_S_K_TO_W_M_K
    k_e = k_e * MCAL_CM_S_K_TO_W_M_K
    # Pr is already unitless

    d_e = d_t * math.sqrt(eps)

    stations = {
        "c": {"Cp": Cp_c, "mu": mu_c, "k": k_c, "Pr": Pr_c, "d": d_c},
        "t": {"Cp": Cp_t, "mu": mu_t, "k": k_t, "Pr": Pr_t, "d": d_t},
        "e": {"Cp": Cp_e, "mu": mu_e, "k": k_e, "Pr": Pr_e, "d": d_e},
    }

    # CEA uses psia, but for SI-based heat transfer correlations use Pa.
    pc_Pa = pc * PSIA_TO_PA

    # NOTE: This assumes SI-consistent inputs: d_* in m, Cp in J/kg-K, mu in Pa·s, pc in Pa, cstar in m/s.
    # If you later add the full Bartz sigma/curvature terms, do it here.
    h_g = {}
    for s, p in stations.items():
        h_g[s] = (
            (0.026 / (d_t ** 0.2))
            * ((p["mu"] ** 0.2) * p["Cp"] / (p["Pr"] ** 0.6))
            * ((pc_Pa / cstar) ** 0.8)
            * ((d_t / p["d"]) ** 1.8)
        )

    # At this point you have:
    # h_g["c"], h_g["t"], h_g["e"]   in ~W/(m^2·K) (with the SI assumption above)
    # and the gas properties per station in `stations`.

    # TODO: continue with wall heating (q"), wall temperature Tw, conduction, stress/thickness, etc.

    return h_g, stations
