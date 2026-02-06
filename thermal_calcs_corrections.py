import math
from rocketcea.cea_obj import CEA_Obj
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from thermal_calcs import *

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

            #for the nozzle, we should section into pieces and calculate flux at each location along. Not neccessary probably
            #for our combustion chamber, as area is constant. Maybe do a reimann sum of small cylinders to approximate area along nozzle and throat
            cp *= 4180  # keeping your conversion factor as-written (verify units!)

            C = rho * cp * delta_x      #flux is through area?
            G = k/delta_x           #thermal diffusivity where?

            h_g = calculate_heat_transfer_coefficient(T_1,T_ad,pc,c_star=c_star,A_t=A_t_m2,A=slice_area,D_t=D_t,r_c=1,M=M,oxname=oxname,fuelname=fuelname,of=of,eps=eps,location=location) #fix all set values 
            q_in = h_g*(T_aw-T_1)
            q_cond = G*(T_1-T_2)    #this is equation for steady state through a plane wall. We need transient through a cylindrical wall, which is a different equation
            q_rad = emissivity * boltzmann * (T_2 ** 4 - T_amb ** 4)    #account for heated gas to wall!
            T_1 = T_1 + (delta_t / C) * (q_in - q_cond)     #q_in - q_cond + q_rad
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