import math

def calculate_wall_temperature(T_initial, delta_t):
    print("Hello World")
    
def calculate_adiabatic_wall_temp(T_ad, gamma, M_e):
    """
    Finds the adiabatic wall temperature given the adiabatic flame temperature,
    the gamma value for the propellants, and the local mach number
    
    :param T_ad: adiabatic flame temperature (K)
    :param gamma: specific heat ratio
    :param M_e: local mach number
    """
    T_e = T_ad / (1 + (gamma - 1)/2 * M_e**2) # Free-stream temperature in Kelvin
    mu_g = 0 # FIND FROM ROCKETCEA
    c_pg = 0 # FIND FROM ROCKETCEA
    k_g = 0 # FIND FROM ROCKETCEA
    Pr = mu_g * c_pg / k_g # Prandtl number of gases 
    r = Pr**(1/3) # Recovery factor
    T_aw = T_e + r(T_ad - T_e) # Adiabatic wall temperature (K)
    return T_aw

def calculate_heat_transfer_coefficient(T_guess, T_aw, T_ad, pc, c_star, A_t, A, D_t, r_c, gamma, M, m = 0.6):
    """
    Finds the gas-side heat transfer coefficient h_g
    
    :param T_guess: The guessed wall temperature (K)
    :param T_aw: Adiabatic wall temperature (K)
    :param T_ad: Adiabatic flame temperature (K)
    :param pc: Chamber pressure (Pa)
    :param c_star: Characteristic velocity (m/s)
    :param A_t: Throat area (m^2)
    :param A: Area at segment (m^2)
    :param D_t: Throat diameter (m)
    :param r_c: Throat radius of curvature (m)
    :param gamma: Ratio of specific heats
    :param M: Local mach
    :param m: Constant
    """
    k_g = 0 # FIND FROM ROCKETCEA
    mu_g = 0 # FIND FROM ROCKETCEA
    c_pg = 0 # FIND FROM ROCKETCEA
    Pr = c_pg * mu_g / k_g # Prandtl number of gases
    sigma = (0.5 * (T_guess/T_ad) * (1 + (gamma - 1)/2 * M**2) + 0.5)**(-(0.8 - m/5)) * (1 + (gamma - 1)/2 * M**2)**(m/5) # Correction factor
    h_g = 0.026/D_t**0.2 * (mu_g**0.2 * c_pg)/(Pr**0.6) * (pc/c_star)**0.8 * (A_t/A)**0.9 * (D_t/r_c)**0.1 * sigma # Convective heat transfer coefficient (W/m^2-K)
    return h_g
