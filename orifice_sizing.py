import numpy as np
import math as math
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matplotlib.pyplot as plt
import throat_sizing

# ----------------------------
# Custom fuel: Paraffin/Al/C40H82, 80/10/10 by mass
# ----------------------------
paraffin_al_c40_card = """
fuel Paraffin80   C 73.0  H 148.0   wt%=80.0
h,cal=-4.4464E+05     t(k)=298.15

fuel Aluminum10   AL 1.0             wt%=10.0
h,cal=0.0             t(k)=298.15

fuel C40H82_10    C 40.0  H 82.0     wt%=10.0
h,cal=-2.0768E+05     t(k)=298.15
"""

add_new_fuel("Paraffin_Al_C40_80_10_10", paraffin_al_c40_card)

# ----------------------------
# Globals / design parameters
# ----------------------------
oxName = 'N2O'
fuelName = 'Paraffin_Al_C40_80_10_10'

pamb = 14.7  # psia (unused here, but keeping)
m_sq_to_in_sq = 1550  # 1 m^2 = 1550 in^2 (unused here)

# Densities (design values)
rhoN2O  = 750.0  # kg/m^3, liquid N2O at ~20–25 C
rhofuel = 950.0  # kg/m^3, paraffin/Al/C40 mixture

# Gas gamma for N2O if needed later (not used in incompressible formula)
g = 1.26  # ox

# Pressures (Pa)
P0 = 750 * 6894.76  # Pa, N2O tank pressure (example)
P1 = P0             # Pa, fuel feed pressure (if same as N2O tank)
P2 = 300 * 6894.76  # Pa, chamber pressure

# Mass flow assumptions
mdot = 0.014279838482639268  # kg/s, oxidizer mass flow (N2O)

m_dot_ox = mdot

# Chamber geometry (inches)
L_star = 45.0  # in
D_c    = 1.25   # in

# Discharge coefficient assumption for orifices
Cd_default = 0.7


def orifice_area(A_t_in2, of, Cd=Cd_default):
    """
    A_t_in2 : throat area in in^2 (used only for chamber geometry / L* calcs)
    Cd      : discharge coefficient (same for fuel and oxidizer here)
    """
    C = CEA_Obj(oxName=oxName, fuelName=fuelName)

    # m dot of fuel
    m_dot_fuel = m_dot_ox / of

    # ----------------------------
    # Chamber geometry from L*
    # ----------------------------
    print("Orifice Sizing Parameters")
    print("----------------------------------")

    print(f"L* for N2O / Paraffin mix: {L_star:.2f} in")

    # A_t is in^2 and L* in in, so V_c is in^3
    V_c = A_t_in2 * L_star
    print(f"Chamber Volume: {V_c:.3f} in^3")
    print(f"Chamber Diameter: {D_c:.3f} in")

    A_c = math.pi * D_c**2 / 4.0
    print(f"Chamber Area: {A_c:.3f} in^2")

    # Add 10% extra length margin (factor 1.1)
    L_c = V_c / (1.1 * A_c)
    print(f"Chamber Length: {L_c:.3f} in")

    # ----------------------------
    # Fuel orifice (incompressible)
    # ----------------------------
    # Pressure drop across fuel orifice
    DeltaP_fuel = P1 - P2  # Pa

    if DeltaP_fuel <= 0:
        raise ValueError("Fuel ΔP must be positive (P1 > P2). Check your pressures.")

    # Cd*A for fuel (m^2)
    CA_fuel = m_dot_fuel / math.sqrt(2.0 * rhofuel * DeltaP_fuel)
    # Actual area (m^2)
    A_fuel = CA_fuel / Cd
    # Convert to mm^2
    A_fuel_mm2 = A_fuel * 1e6

    print(f"Fuel mass flow (Paraffin mix): {m_dot_fuel:.6f} kg/s")
    print(f"Fuel Orifice C_d*A: {CA_fuel:.6e} m^2")
    print(f"Fuel Orifice Area (A, C_d={Cd}): {A_fuel_mm2:.3f} mm^2")

    # ----------------------------
    # N2O orifice (incompressible)
    # ----------------------------
    # Pressure drop across N2O orifice
    DeltaP_ox = P0 - P2  # Pa

    if DeltaP_ox <= 0:
        raise ValueError("N2O ΔP must be positive (P0 > P2). Check your pressures.")

    # Cd*A for N2O (m^2)
    CA_ox = m_dot_ox / math.sqrt(2.0 * rhoN2O * DeltaP_ox)
    # Actual area (m^2)
    A_ox = CA_ox / Cd
    # Convert to mm^2
    A_ox_mm2 = A_ox * 1e6

    print(f"N2O mass flow: {m_dot_ox:.6f} kg/s")
    print(f"N2O Orifice C_d*A: {CA_ox:.6e} m^2")
    print(f"N2O Orifice Area (A, C_d={Cd}): {A_ox_mm2:.3f} mm^2")

    print("----------------------------------")

    # Return areas in case you want to programmatically use them
    return {
        "A_fuel_m2": A_fuel,
        "A_fuel_mm2": A_fuel_mm2,
        "A_ox_m2": A_ox,
        "A_ox_mm2": A_ox_mm2,
        "L_c_in": L_c,
        "V_c_in3": V_c
    }