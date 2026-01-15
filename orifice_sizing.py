import math as math
import traceback

def orifice_area(A_t_in2, of, Cd, mdot, L_star, D_c, D_t, D_e, conv_angle, div_angle, P1, P2, rho_ox):
    #print("\n[DEBUG] orifice_area called from:")
    #traceback.print_stack(limit=4)

    """
    Finds the sizes of various orifices for the engine, include the size of the engine itself
    
    :param A_t_in2: Throat area in in^2
    :param of: of ratio
    :param Cd: fuel coefficient
    :param mdot: mass flow rate (kg/s)
    :param L_star: characteristic length (in)
    :param D_c: chamber diameter (in)
    :param D_t: throat diameter (in)
    :param D_e: exit diameter (in)
    :param conv_angle: convergence angle (degrees)
    :param div_angle: divergence angle (degrees)
    :param P1: inlet pressure (pa)
    :param P2: chamber pressure (pa)
    :param rho_ox: density of oxidizer (kg/m^3)
    """

    # m dot of N2O
    m_dot_ox = mdot * (of / (1 + of))

    # ----------------------------
    # Chamber geometry from L*
    # ----------------------------
    print("----------------------------------")
    print("Orifice Sizing Parameters")
    print("----------------------------------")

    print(f"L* for N2O / Paraffin mix: {L_star:.2f} in")

    # A_t is in^2 and L* in in, so V_c is in^3
    V_c = A_t_in2 * L_star
    print(f"Chamber Volume: {V_c:.3f} in^3")
    print(f"Chamber Diameter: {D_c:.3f} in")

    A_c = math.pi * D_c**2 / 4.0
    print(f"Chamber Area: {A_c:.3f} in^2")

    # Treat converging section as a cone with 30 degree angle
    L_conv = (D_c - D_t) / (2 * math.tan(math.radians(conv_angle)))
    V_conv = math.pi * L_conv / 3 * ((D_c/2)**2 + (D_c * D_t / 4) + (D_t/2)**2)

    # Find the length of the straight section from this
    L_straight = (V_c - V_conv)/ A_c
    
    print(f"Chamber Length (Straight): {L_straight:.3f} in")
    print(f"Convergence Angle: {conv_angle} degrees")
    print(f"Chamber Length (Converging): {L_conv:.3f} in")

    # Length of diverging section for conical nozzle
    L_div = (D_e - D_t) / (2 * math.tan(math.radians(div_angle)))  
    print(f"Divergence Angle: {div_angle} degrees")
    print(f"Length of Diverging Section (Straight Nozzle: {L_div:.3f} in")

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

    return(A_ox, A_ox_mm2, L_straight, L_conv, V_c, L_div)
   