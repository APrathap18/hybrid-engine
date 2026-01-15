import math
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PA_TO_PSI = 0.000145038
cal_to_kj = 4184

def calculate_wall_temperature(T_initial, delta_t,T_ad,t_wall,pc,oxname,fuelname,of,eps,mach_exit,burn_time,D_c,c_star,D_t,D_e,emissivity = 0.3, boltzmann = 5.67e-8, T_amb = 298):
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
    
    T_aw_list = []
    locations = ['c', 't', 'e']
    #X_dict = {i:  round((i*0.1),1) for i in range(0, int(l / 0.1) + 2)}

    T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,0,oxname,fuelname,of,eps,location='c'))
    T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,1,oxname,fuelname,of,eps,location='t'))
    T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,mach_exit,oxname,fuelname,of,eps,location='e'))

    # for X in X_dict:
    #     if X < l_c:
    #         M_e = 0 #nozzle geometry
    #         T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='c'))
    #     elif X >= l_c and X < l:
    #         M_e = 1 #nozzle geometry
    #         T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='t'))
    #     elif X == l:
    #         M_e = mach_exit #nozzle geometry
    #         T_aw_list.append(calculate_adiabatic_wall_temp(T_ad,pc,M_e,oxname,fuelname,of,eps,location='e'))
    N = 2
    delta_x = t_wall/N

    # material property tables (keys are K rounded to nearest 100)
    rho_dict = {300:7894,400:7860,500:7823,600:7783,700:7742,800:7698,900:7652,1000:7603,1100:7552,1200:7499,1300:7444,1400:7386,1500:7326,1600:7264,1700:7199}
    c_p_dict = {300:0.1219,400:0.1251,500:0.1283,600:0.1315,700:0.1348,800:0.1380,900:0.1412,1000:0.1444,1100:0.1476,1200:0.1509,1300:0.1541,1400:0.1573,1500:0.1605,1600:0.1638,1700:0.1670,1800:0.1900}
    k_dict   = {300:12.97,400:14.59,500:16.20,600:17.82,700:19.44,800:21.06,900:22.67,1000:24.29,1100:25.91,1200:27.53,1300:29.14,1400:30.76,1500:32.38,1600:34.00,1700:35.61,1800:18.14}  # fixed obvious typos
    
    A_t_m2 = math.pi * (D_t / 2)**2

    time_array = np.arange(0, burn_time + 1e-12, delta_t)

    chamber_wall_hot_temp = []
    throat_wall_hot_temp = []
    exit_wall_hot_temp = []

    chamber_wall_outer_temp = []
    throat_wall_outer_temp = []
    exit_wall_outer_temp = []

    for location in locations:
        T_1 = T_initial
        T_2 = T_initial
        T_aw = T_aw_list[locations.index(location)]
        
            
        if location == 'c':
            slice_area = math.pi * (D_c / 2)**2
            M = 0
        elif location == 't':
            slice_area = math.pi * (D_t / 2)**2
            M = 1
        elif location == 'e':
            slice_area = math.pi * (D_e / 2)**2
            M = mach_exit
            
        for t in time_array:
            n_T = int(round(T_1 / 100.0) * 100) 
            n_T = max(300, min(n_T, 1600)) 
            n_T2 = n_T + 100

            frac = (T_1 - n_T) / 100.0  # 0..1 inside the bin

            k   = k_dict[n_T]   + frac * (k_dict[n_T2]   - k_dict[n_T])
            rho = rho_dict[n_T] + frac * (rho_dict[n_T2] - rho_dict[n_T])
            cp  = c_p_dict[n_T] + frac * (c_p_dict[n_T2] - c_p_dict[n_T])

            cp *= 4180  # keeping your conversion factor as-written (verify units!)

            C = rho * cp * delta_x
            G = k/delta_x

            h_g = calculate_heat_transfer_coefficient(T_1,T_ad,pc,c_star=c_star,A_t=A_t_m2,A=slice_area,D_t=D_t,r_c=1,M=M,oxname=oxname,fuelname=fuelname,of=of,eps=eps,location=location) #fix all set values 
            q_in = h_g*(T_aw-T_1)
            q_cond = G*(T_1-T_2)
            q_rad = emissivity * boltzmann * (T_2 ** 4 - T_amb ** 4)
            T_1 = T_1 + (delta_t / C) * (q_in - q_cond)
            T_2 = T_2 + (delta_t / C) * (q_cond - q_rad)

            if location == 'c':
                chamber_wall_hot_temp.append(T_1)
                chamber_wall_outer_temp.append(T_2)
            elif location == 't':
                throat_wall_hot_temp.append(T_1)
                throat_wall_outer_temp.append(T_2)
            elif location == 'e':
                exit_wall_hot_temp.append(T_1)
                exit_wall_outer_temp.append(T_2)
    df = pd.DataFrame({
        "Time (s)": time_array,
        "Chamber Hot-Wall Temperature (K)": chamber_wall_hot_temp,
        "Throat Hot-Wall Temperature (K)": throat_wall_hot_temp,
        "Exit Hot-Wall Temperature (K)": exit_wall_hot_temp,
    })

    print(df.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    plot_thermal_data(time_array, chamber_wall_hot_temp, throat_wall_hot_temp, exit_wall_hot_temp, chamber_wall_outer_temp, throat_wall_outer_temp, exit_wall_outer_temp, oxname, fuelname, pc, eps)

    # for T_w in T_array:
    #     n_T = int(round(T_w / 100.0) * 100)
    #     n_T = max(300, min(n_T, 1700))
    #     n_T2 = n_T + 100

    #     frac = (T_w - n_T) / 100.0  # 0..1 inside the bin

    #     k   = k_dict[n_T]   + frac * (k_dict[n_T2]   - k_dict[n_T])
    #     rho = rho_dict[n_T] + frac * (rho_dict[n_T2] - rho_dict[n_T])
    #     cp  = c_p_dict[n_T] + frac * (c_p_dict[n_T2] - c_p_dict[n_T])
    #     cp *= 4180  # keeping your conversion factor as-written (verify units!)

    #     print(f"Guessed Wall Temperature: {T_w} K")

    #     time_array = np.arange(0, burn_time + 1e-12, delta_t)
        
    #     for location in locations:
    #         chamber_wall_temp = []
    #         throat_wall_temp = []
    #         exit_wall_temp = []

    #         T_1 = T_initial
    #         T_2 = T_initial
    #         T_aw = T_aw_list[locations.index(location)]
    #         C = rho * cp * delta_x
    #         G = k/delta_x
            
    #         if location == 'c':
    #             slice_area = math.pi * (D_c / 2)**2
    #             M = 0
    #         elif location == 't':
    #             slice_area = math.pi * (D_t / 2)**2
    #             M = 1
    #         elif location == 'e':
    #             slice_area = math.pi * (D_e / 2)**2
    #             M = mach_exit
            
    #         for t in time_array:
    #             h_g = calculate_heat_transfer_coefficient(T_1,T_ad,pc,c_star=c_star,A_t=A_t_m2,A=slice_area,D_t=D_t,r_c=1,M=M,oxname=oxname,fuelname=fuelname,of=of,eps=eps,location=location) #fix all set values 
    #             q_in = h_g*(T_aw-T_1)
    #             q_cond = G*(T_1-T_2)
    #             T_1 = T_1 + (delta_t / C) * (q_in - q_cond)
    #             T_2 = T_2 + (delta_t / C) * q_cond

    #             if location == 'c':
    #                 chamber_wall_temp.append(T_1)
    #             elif location == 't':
    #                 throat_wall_temp.append(T_1)
    #             elif location == 'e':
    #                 exit_wall_temp.append(T_1)
            
    #     df = pd.DataFrame({
    #         "Time (s)": time_array,
    #         "Chamber Temperature (K)": chamber_wall_temp,
    #         "Throat Temperature (K)": throat_wall_temp,
    #         "Exit Temperature (K)": exit_wall_temp,
    #     })

    #     print(df.to_string(index=False, float_format=lambda x: f"{x:.1f}"))

        
        
def plot_thermal_data(time_array, chamber_temp, throat_temp, exit_temp, chamber_outer, throat_outer, exit_outer, oxName, fuelName, pc, eps):
    # creates a 3 x 3 subplot
    fig, axes = plt.subplots(3, 3, figsize=(10, 8), sharex=True)

    # flatten axes
    axes = axes.ravel()

    # graph title
    fig.suptitle(f'{oxName}/{fuelName} @ Pc = {pc:.1f} psia, eps = {eps:.1f}', fontsize=14, fontweight="bold")

    # suplot 1, chamber temp inner wall
    axes[0].plot(time_array, chamber_temp, marker='o', linestyle='-')
    axes[0].set_ylabel('Chamber Inner Wall Temp [K]')
    axes[0].grid(True)
    axes[0].set_title('Chamber Inner Wall Temp [K] vs. Time [s]')

    # suplot 2, chamber temp outer wall
    axes[1].plot(time_array, chamber_outer, marker='o', linestyle='-')
    axes[1].set_ylabel('Chamber Outer Wall Temp [K]')
    axes[1].grid(True)
    axes[1].set_title('Chamber Outer Wall Temp [K] vs. Time [s]')

    chamber_delta_T = [a - b for a, b in zip(chamber_temp, chamber_outer)]

    # suplot 3, chamber delta T
    axes[2].plot(time_array, (chamber_delta_T), marker='o', linestyle='-')
    axes[2].set_ylabel('Temp Difference [K]')
    axes[2].grid(True)
    axes[2].set_title('Chamber Inner-Outer Wall Temp Difference [K] vs. Time [s]')

    # suplot 4, throat temp inner wall
    axes[3].plot(time_array, throat_temp, marker='o', linestyle='-')
    axes[3].set_ylabel('Throat Inner Wall Temp [K]')
    axes[3].grid(True)
    axes[3].set_title('Throat Inner Wall Temp [K] vs. Time [s]')

    # suplot 5, throat temp outer wall
    axes[4].plot(time_array, throat_outer, marker='o', linestyle='-')
    axes[4].set_ylabel('Throat Outer Wall Temp [K]')
    axes[4].grid(True)
    axes[4].set_title('Throat Outer Wall Temp [K] vs. Time [s]')

    throat_delta_T = [a - b for a, b in zip(throat_temp, throat_outer)]

    # suplot 6, throat delta T
    axes[5].plot(time_array, throat_delta_T, marker='o', linestyle='-')
    axes[5].set_ylabel('Temp Difference [K]')
    axes[5].grid(True)
    axes[5].set_title('Throat Inner-Outer Wall Temp Difference [K] vs. Time [s]')

    # suplot 7, exit temp inner wall
    axes[6].plot(time_array, exit_temp, marker='o', linestyle='-')
    axes[6].set_ylabel('Exit Inner Wall Temp [K]')
    axes[6].grid(True)
    axes[6].set_title('Exit Inner Wall Temp [K] vs. Time [s]')

    # suplot 8, exit temp outer wall
    axes[7].plot(time_array, exit_outer, marker='o', linestyle='-')
    axes[7].set_ylabel('Exit Outer Wall Temp [K]')
    axes[7].grid(True)
    axes[7].set_title('Exit Outer Wall Temp [K] vs. Time [s]')

    exit_delta_T = [a - b for a, b in zip(exit_temp, exit_outer)]

    # suplot 9, exit delta T
    axes[8].plot(time_array, exit_delta_T, marker='o', linestyle='-')
    axes[8].set_ylabel('Temp Difference [K]')
    axes[8].grid(True)
    axes[8].set_title('Exit Inner-Outer Wall Temp Difference [K] vs. Time [s]')

    # show plot
    plt.tight_layout()
    plt.show()


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
    pc_psi = pc * PA_TO_PSI

    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)
    _, gamma = cea_obj.get_Chamber_MolWt_gamma(pc_psi, of, eps)

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

def calculate_heat_transfer_coefficient(T_guess, T_ad, pc, c_star, A_t, A, D_t, r_c, M, oxname, fuelname, of, eps, location, m = 0.6):
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
    pc_psi = pc * PA_TO_PSI
    cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)
    _, gamma = cea_obj.get_Chamber_MolWt_gamma(pc_psi, of, eps)

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
