###############################
# Ablative Calcs Mk. 1
###############################



import math
import numpy as np
from thermal_calcs import *
from rocketcea.cea_obj import CEA_Obj



def ablative_rate(T_initial, delta_t,T_ad,T_0g,t_wall,pc,oxname,fuelname,of,eps,mach_exit,burn_time,D_c,c_star,D_t,D_e,emissivity = 0.3, boltzmann = 5.67e-8, T_amb = 298, abl_t0 = 5):
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
    rho_abl_dict = {}
    cp_abl_dict = {}
    k_abl_dict = {}
    rec_abl_dict = {}

    A_t_m2 = math.pi * (D_t / 2)**2

    time_array = np.arange(0, burn_time + 1e-12, delta_t)

    chamber_wall_hot_temp = []
    throat_wall_hot_temp = []
    exit_wall_hot_temp = []

    chamber_wall_outer_temp = []
    throat_wall_outer_temp = []
    exit_wall_outer_temp = []

    chamber_abl_temp = []
    throat_abl_temp =[]
    exit_abl_temp = []

    for location in locations:
        T_char = T_initial
        T_v = T_initial
        T_1 = T_initial
        T_2 = T_initial
        

        T_aw = T_aw_list[locations.index(location)]
        abl_t = abl_t0
        char_t = 0
        if location == 'c':
            slice_area = math.pi * (D_c / 2)**2
            M = 0
            ablated = True
            charred = True
        elif location == 't':
            slice_area = math.pi * (D_t / 2 + abl_t)**2
            M = 1
            ablated = False
            charred = False
        elif location == 'e':
            slice_area = math.pi * (D_e / 2 + abl_t)**2
            M = mach_exit
            ablated = False
            charred = False
            
        for t in time_array:
            n_T = int(round(T_1 / 100.0) * 100) 
            n_T = max(300, min(n_T, 1600)) 
            n_T2 = n_T + 100
        
            

            frac = (T_1 - n_T) / 100.0  # 0..1 inside the bin

            k_bod   = k_dict[n_T]   + frac * (k_dict[n_T2]   - k_dict[n_T])
            rho_bod = rho_dict[n_T] + frac * (rho_dict[n_T2] - rho_dict[n_T])
            cp_bod  = c_p_dict[n_T] + frac * (c_p_dict[n_T2] - c_p_dict[n_T])

            cp_bod *= 4180  # keeping your conversion factor as-written (verify units!)

            C = rho_bod * cp_bod * delta_x
            G = k_bod/delta_x
            C_char = rho_char*cp_char*char_t
            
            if charred == False:            
                h_g = calculate_heat_transfer_coefficient(T_1,T_ad,pc,c_star=c_star,A_t=A_t_m2,A=slice_area,D_t=D_t,r_c=1,M=M,oxname=oxname,fuelname=fuelname,of=of,eps=eps,location=location)
                q_conv = h_g*(T_aw - T_v)
                m_rate = m_loss_rate_dict[int(round(q_conv,-2))]
                abl_rate = m_rate/rho_abl
                abl_t -= abl_rate*delta_t
                char_t = abl_t0 - abl_t
                T_char = max(0, (T_char + delta_t/(rho_char*cp_char*char_t)*(q_conv-(m_rate*H_abl))))
                charred == True




            if ablated == False:
                h_g = calculate_heat_transfer_coefficient(T_1,T_ad,pc,c_star=c_star,A_t=A_t_m2,A=slice_area,D_t=D_t,r_c=1,M=M,oxname=oxname,fuelname=fuelname,of=of,eps=eps,location=location) #fix all set values 
                q_in = h_g*(T_aw-T_1)
                q_cond1 = G*(T_1-T_2)
                q_cond2 = G*(T_2-T_3)
                q_rad = emissivity * boltzmann * (T_ad**4 - T_1**4)
                q_rad2 = emissivity * boltzmann * (T_2 ** 4 - T_amb ** 4)
                T_1 = T_1 + (delta_t / C) * (q_in - q_cond1 + q_rad) - h_abl/cp_abl
                T_2 = T_2 + (delta_t / C) * (q_cond1 - q_cond2)
                T_3 = T_3 + (delta_t / C) * (q_cond2)
        

            if location == 'c':
                chamber_wall_hot_temp.append(T_1)
                chamber_wall_outer_temp.append(T_2)
                chamber_abl_temp.append('none')
            elif location == 't' and ablated == True:
                throat_abl_temp.append('none')
                throat_wall_hot_temp.append(T_1)
                throat_wall_outer_temp.append(T_2)
                slice_area = math.pi * (D_t / 2)**2
            elif location == 't' and ablated == False:
                throat_abl_temp.append(T_1)
                throat_wall_hot_temp.append(T_2)
                throat_wall_outer_temp.append(T_3)
                slice_area = math.pi * (D_t / 2 + abl_t)**2
            elif location == 'e' and ablated == True:
                exit_abl_temp.append('none')
                exit_wall_hot_temp.append(T_1)
                exit_wall_outer_temp.append(T_2)
                slice_area = math.pi * (D_e / 2)**2
            elif location == 'e' and ablated == False:
                exit_abl_temp.append(T_1)
                exit_wall_hot_temp.append(T_2)
                exit_wall_outer_temp.append(T_3)
                slice_area = math.pi * (D_e / 2 + abl_t)**2
            if abl_t <= 0 and ablated == False:
                abl_fail = t
                ablated = True
                T_1 = T_2
                T_2 = T_3
            

            

    df = pd.DataFrame({
        "Time (s)": time_array,
        "Chamber Hot-Wall Temperature (K)": chamber_wall_hot_temp,
        "Throat Hot-Wall Temperature (K)": throat_wall_hot_temp,
        "Exit Hot-Wall Temperature (K)": exit_wall_hot_temp,
    })

    print(df.to_string(index=False, float_format=lambda x: f"{x:.1f}"))
    return abl_fail

