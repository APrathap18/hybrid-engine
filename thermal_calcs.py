import math
from rocketcea.cea_obj import CEA_Obj

PA_TO_PSI = 0.000145038

def calculate_wall_temperature(T_initial, delta_t):
    print("Hello World")
    
def calculate_adiabatic_wall_temp(T_ad, pc, M_e, oxname, fuelname, of, eps, location):
    """
    Finds the adiabatic wall temperature given the adiabatic flame temperature,
    the gamma value for the propellants, and the local mach number
    
    :param T_ad: adiabatic flame temperature (K)
    :param pc: chamber pressure (Pa)
    :param M_e: local mach number
    :param oxname: oxidizer name (string)
    :param fuelname: fuel name (string)
    :param of: O/F ratio
    :param eps: Expansion ratio
    :param location: String indicating if the function should be solved for the chamber ("c"), throat ("t"), or exit ("e")
    """
    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)

    pc_psi = pc * PA_TO_PSI

    T_e = T_ad / (1 + (gamma - 1)/2 * M_e**2) # Free-stream temperature in Kelvin

    if location == 'c':
        _, gamma = cea_obj.get_Chamber_MolWt_gamma(Pc = pc_psi, MR = of, eps = eps)
        _, _, _, Pr = cea_obj.get_Chamber_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    elif location == 't':
        _, _, _, Pr = cea_obj.get_Throat_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    elif location == 'e':
        _, _, _, Pr = cea_obj.get_Exit_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    else:
        raise ValueError("Invalid location provided. Specify if values should be calculated at chamber, throat, or exit")
    
    r = Pr**(1/3) # Recovery factor
    T_aw = T_e + r * (T_ad - T_e) # Adiabatic wall temperature (K)
    return T_aw

def calculate_heat_transfer_coefficient(T_guess, T_ad, pc, c_star, A_t, A, D_t, r_c, gamma, M, oxname, fuelname, of, eps, location, m = 0.6):
    """
    Finds the gas-side heat transfer coefficient h_g
    
    :param T_guess: The guessed wall temperature (K)
    :param T_ad: Adiabatic flame temperature (K)
    :param pc: Chamber pressure (Pa)
    :param c_star: Characteristic velocity (m/s)
    :param A_t: Throat area (m^2)
    :param A: Area at segment (m^2)
    :param D_t: Throat diameter (m)
    :param r_c: Throat radius of curvature (m)
    :param gamma: Ratio of specific heats
    :param M: Local mach
    :param oxname: oxidizer name (string)
    :param fuelname: fuel name (string)
    :param of: O/F ratio
    :param eps: Expansion ratio
    :param location: String indicating if the function should be solved for the chamber ("c"), throat ("t"), or exit ("e")
    :param m: Constant
    """
    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)

    pc_psi = pc * PA_TO_PSI

    if location == 'c':
        cp, mu, _, Pr = cea_obj.get_Chamber_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    elif location == 't':
        cp, mu, _, Pr = cea_obj.get_Throat_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    elif location == 'e':
        cp, mu, _, Pr = cea_obj.get_Exit_Transport(Pc = pc_psi, MR = of, eps = eps, frozen = 0)
    else:
        raise ValueError("Invalid location provided. Specify if values should be calculated at chamber, throat, or exit")

    # Cp: cal/(g·K) -> J/(kg·K)
    cp_SI = cp * 4184.0

    # viscosity: millipoise -> Pa·s
    # 1 P = 0.1 Pa·s, 1 mP = 1e-3 P = 1e-4 Pa·s
    mu_SI = mu * 1e-4

    sigma = (0.5 * (T_guess/T_ad) * (1 + (gamma - 1)/2 * M**2) + 0.5)**(-(0.8 - m/5)) * (1 + (gamma - 1)/2 * M**2)**(m/5) # Correction factor
    h_g = 0.026/D_t**0.2 * (mu_SI**0.2 * cp_SI)/(Pr**0.6) * (pc/c_star)**0.8 * (A_t/A)**0.9 * (D_t/r_c)**0.1 * sigma # Convective heat transfer coefficient (W/m^2-K)
    return h_g
