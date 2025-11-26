import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matplotlib.pyplot as plt


# ----------------------------
# Custom fuel: Paraffin/Al/C40H82, 80/10/10 by mass
# ----------------------------
paraffin_al_c40_card = """
fuel Paraffin80   C 73.0  H 148.0   wt%=80.0
h,cal=-4.4464E+05     t(k)=298.15

fuel Aluminum10  AL 1.0             wt%=10.0
h,cal=0.0             t(k)=298.15

fuel C40H82_10   C 40.0  H 82.0     wt%=10.0
h,cal=-2.0768E+05     t(k)=298.15
"""

add_new_fuel("Paraffin_Al_C40_80_10_10", paraffin_al_c40_card)

# ----------------------------
# Globals / design parameters
# ----------------------------
oxName = 'N2O'
fuelName = 'Paraffin_Al_C40_80_10_10'
pamb = 14.7 # psia
m_sq_to_in_sq = 1550 #1 m^2 = 1550 in^2
n2o_density = 750 # kg/m^3
fuel_density = 950 # kg/m^3
eps = 4

# This is a function that takes O/F ratio, chamber pressure (psia), 
# and thrust (N) as inputs, and outputs the diameter (in) of the throat 
def throat_sizing_function(of, pc, F):

    C = CEA_Obj(oxName=oxName, fuelName=fuelName)

    # 1) Get vacuum Isp and C* from CEA (quiet)
    Isp_vac, Cstar_ft_s, Tc, MW, gamma = C.get_IvacCstrTc_ChmMwGam(
        Pc=pc,
        MR=of,
        eps=eps       # your chosen area ratio (Ae/At)
    )

    # 2) Get "CEA Isp" at this eps (also quiet)
    Isp_CEA = C.get_Isp(Pc=pc, MR=of, eps=eps)  # [s]

    # 3) Ambient correction (no CalcPCoPE used)
    Pamb   = 14.7      # psia
    g0_ft  = 32.174    # ft/s^2

    Isp_amb = Isp_CEA - Cstar_ft_s * Pamb * eps / (pc * g0_ft)  # [s]

    # 4) Convert to exhaust velocity [m/s]
    v_e_ft_s = Isp_amb * g0_ft       # ft/s
    v_e      = v_e_ft_s * 0.3048     # m/s

    # Optional: assume some loss vs ideal if you want, e.g. 90%
    v_e_eff = 0.9 * v_e

    # 5) Total mass flow from thrust
    mdot_total = F / v_e_eff         # kg/s

    # Split into ox/fuel using O/F
    m_dot_ox   = mdot_total * of / (1.0 + of)
    m_dot_fuel = mdot_total / (1.0 + of)

    # 6) Compute throat area from Pc and C*
    pc_pa = pc * 6894.76             # psia -> Pa
    Cstar_m_s = Cstar_ft_s * 0.3048  # ft/s -> m/s

    At = mdot_total * Cstar_m_s / pc_pa   # m^2

    Dt_m  = 2.0 * np.sqrt(At / np.pi)
    Dt_in = Dt_m / 0.0254                 # 1 in = 0.0254 m

    print("----------------------------------")
    print("Throat Sizing Parameters")
    print("----------------------------------")
    print(f"gamma                : {gamma:.4f}")
    print(f"Isp_vac (s)         : {Isp_vac:.2f}")
    print(f"Isp_CEA (s)         : {Isp_CEA:.2f}")
    print(f"Isp_amb (s)         : {Isp_amb:.2f}")
    print(f"Exhaust velocity (m/s): {v_e_eff:.2f}")
    print(f"Total m_dot (kg/s)  : {mdot_total:.4f}")
    print(f"N2O m_dot (kg/s)    : {m_dot_ox:.4f}")
    print(f"Fuel m_dot (kg/s)   : {m_dot_fuel:.4f}")
    print(f"Throat Area (m^2)   : {At:.6e}")
    print(f"Throat Area (in^2)  : {At * m_sq_to_in_sq:.4f}")
    print(f"Throat Diameter (in): {Dt_in:.3f}")
    print("----------------------------------")

    return(At * m_sq_to_in_sq)