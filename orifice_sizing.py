import math as math
import traceback

def orifice_area(A_t_in2, of, Cd, mdot, L_star, D_c, P1, P2, rho_ox, rho_fuel, pamb, oxName, fuelName):
    #print("\n[DEBUG] orifice_area called from:")
    #traceback.print_stack(limit=4)

    """
    A_t_in2 : throat area in in^2 (used only for chamber geometry / L* calcs)
    Cd      : discharge coefficient (same for fuel and oxidizer here)
    """

    # m dot of N2O
    m_dot_ox   = mdot * (of / (1 + of))

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
    # Convert to in^2
    A_ox_in2 = A_ox_mm2 * 0.00155

    print(f"N2O Orifice Area (A, C_d={Cd}): {A_ox_mm2:.3f} mm^2")
    print(f"N2O Orifice Area (A, C_d={Cd}): {A_ox_in2:.3f} in^2")

    print("----------------------------------")

    # Return areas in case you want to programmatically use them
    return {
        "A_ox_m2": A_ox,
        "A_ox_mm2": A_ox_mm2,
        "L_c_in": L_c,
        "V_c_in3": V_c
    }