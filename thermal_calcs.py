import math
from rocketcea.cea_obj import CEA_Obj
import numpy as np
PA_TO_PSI = 0.000145038
cal_to_kj = 4184

def calculate_wall_temperature(T_initial, delta_t,T_ad,t,l,l_c,pc,oxname,fuelname,of,eps):
    """
    Calculates the wall temperature of the combustion chamber
    
    :param T_initial: Initial wall temperature (K)
    :param delta_t: Time step (s)
    :param T_ad: Adiabatic flame temperature (K)
    :param t: Thickness (m)
    :param l: Total length (m)
    :param l_c: Chamber length (m)
    :param pc: Chamber pressure (Pa)
    :param oxname: Oxidizer name (string)
    :param fuelname: Fuel type (string)
    :param of: O/F ratio
    :param eps: Expansion ratio
    """
    
    T_array = np.arrange(T_initial,T_ad,100)
    T_aw_list = []
    X_dict = {i:  round((i*0.1),1) for i in range(0, int(l / 0.1) + 2)}
    for X in X_dict:
        if X_dict(X) < l_c:
            M_e = 0 #nozzle geometry
            T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='c'))
        elif X_dict(X) >= l_c and X < l:
            M_e = 1 #nozzle geometry
            T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='t'))
        elif X_dict(X) == l:
            M_e = 3 #nozzle geometry
            T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='e'))
    N = 10
    delta_x = t/N
    rho_dict = {300:7894,400:7860,500:7823,600:7783,700:7742,800:7698,900:7652,1000:7603,1100:7552,1200:7499,1300:7444,1400:7386,1500:7326,1600:7264,1700:7199} #replace with real densities. need to think about melted steel
    c_p_dict = {300:0.1219,400:0.1251,500:0.1283,600:0.1315,700:0.1348,800:0.1380,900:0.1412,1000:0.1444,1100:0.1476,1200:0.1509,1300:0.1541,1400:0.1573,1500:0.1605,1600:0.1638,1700:0.1670,1800:0.1900}
    k_dict = {300:12.97,400:14.59,500:16.2,600:17.82,700:19.44,800:21.06,900:22.67,1000:24.29,1100:25.91,1200:27.53,1300:29.14,1400:30.76,1500:32.38,1600:34.00,1700:3561,1800:1814}
    for T_w in T_array:
        n_T = round(T_w,-2)
        k = k_dict(n_T) + (100)/(k_dict(round((n_T+100),-2)) - k_dict(n_T))*(T_aw - n_T) #assumes roughly linear between steps
        rho = rho_dict(n_T) + (100)/(rho_dict(round((n_T+100),-2)) - rho_dict(n_T))*(T_aw - n_T) #only works for temps under 1700K, where 304SS melts
        c_p = (c_p_dict(n_T) + (100)/(c_p_dict(round((n_T+100),-2)) - c_p_dict(n_T))*(T_aw - n_T))*4180
        for X in X_dict:
            T_1 = T_initial
            T_2 = T_initial
            time_array = np.arrange(0,5,delta_t) #replace 5 with burn_time
            for t in time_array:
                T_aw = T_aw_list(X)
                C = rho * c_p * delta_x
                G = k/delta_x
                h_g = calculate_heat_transfer_coefficient(T_w,T_ad,pc,c_star=1,A_t=1,A=1,D_t=1,r_c=1,gamma=1,M=1,oxname='s',fuelname='s',of='s',eps='s',location='tbd') #fix all set values 
                q_in = h_g*(T_aw-T_w)
                q_cond = G*(T_1-T_2)
                T_1 = T_1 + delta_t/(C*(q_in-q_cond))
                T_2 = T_2 + delta_t/(C*q_cond)
                print(T_1)


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
