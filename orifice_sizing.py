import numpy as np
import math as math
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matplotlib.pyplot as plt
import throat_sizing

def orifice_area(A_t_in2, of, Cd, mdot, L_star, D_c, P1, P2, rho_ox, rho_fuel, pamb, oxName, fuelName):
    """
    A_t_in2 : throat area in in^2 (used only for chamber geometry / L* calcs)
    Cd      : discharge coefficient (same for fuel and oxidizer here)
    """
    C = CEA_Obj(oxName=oxName, fuelName=fuelName)

    # m dot of fuel

    m_dot_fuel = mdot / of
    m_dot_ox = mdot * of

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
    CA_fuel = m_dot_fuel / math.sqrt(2.0 * rho_fuel * DeltaP_fuel)
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
    DeltaP_ox = P1 - P2  # Pa

    if DeltaP_ox <= 0:
        raise ValueError("N2O ΔP must be positive (P0 > P2). Check your pressures.")

    # Cd*A for N2O (m^2)
    CA_ox = m_dot_ox / math.sqrt(2.0 * rho_ox * DeltaP_ox)
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