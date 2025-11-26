import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matplotlib.pyplot as plt

# ----------------------------
# Custom fuel: Paraffin/Al/C40H82, 80/10/10 by mass
# ----------------------------
paraffin_al_c40_card = """
fuel Paraffin80   C 73.0  H 148.0   wt%=80.0
h,cal=-4.4464E+05     t(k)=298.15

fuel Aluminum10  AL 1.0             wt%=10.0
h,cal=0.0             t(k)=298.15

fuel C40H82_10   C 40.0  H 82.0     wt%=10.0
h,cal=-2.0768E+05     t(k)=298.15
"""

add_new_fuel("Paraffin_Al_C40_80_10_10", paraffin_al_c40_card)

# ----------------------------
# Globals / design parameters
# ----------------------------
oxName = 'N2O'
fuelName = 'Paraffin_Al_C40_80_10_10'
pamb = 14.7 # psia
m_sq_to_in_sq = 1550 #1 m^2 = 1550 in^2
n2o_density = 750 # kg/m^3
fuel_density = 950 # kg/m^3

def plot_OF(pc, eps):
    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxName, fuelName = fuelName)

    # lists for chamber temp, throat temp, exhaust temp, isp, and o/f ratios
    comb_temp = []
    throat_temp = []
    exhaust_temp = []
    isp_list = []
    of_list = []

    g0_ft = 32.174  # ft/s^2

    # iterates through o/f ratios in 0.5 step sizes
    for of in np.arange(0.5, 10.0, 0.25):
        try:
            # add o/f ratio to list
            of_list.append(of)

            # get all temperatures in one list (degrees R)
            temp = cea_obj.get_Temperatures(Pc = pc, MR = of, eps = eps)

            # append temps to corresponding lists, converting to Kelvin
            comb_temp.append(temp[0] * 5/9)
            throat_temp.append(temp[1] * 5/9)
            exhaust_temp.append(temp[2] * 5/9)

            # Quiet CEA performance call
            Isp_vac, Cstar_ft_s, Tc, MW, gamma = \
                cea_obj.get_IvacCstrTc_ChmMwGam(Pc=pc, MR=of, eps=eps)

            # Ambient Isp at pamb using RocketCEA's documented correction
            Isp_amb = Isp_vac - Cstar_ft_s * pamb * eps / (pc * g0_ft)

            # Store ambient Isp in seconds
            isp_list.append(Isp_amb)
        except Exception as e:
            # prevents errors
            print(f"Skipping MR = {of:.2f}: {e}")
    
    # creates a 2 x 2 subplot
    fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True)

    # flatten axes
    axes = axes.ravel()

    # graph title
    fig.suptitle(f'{oxName}/{fuelName} @ Pc = {pc} psia, eps = {eps}, Pamb = {pamb} psia', fontsize=14, fontweight="bold")

    # suplot 1, chamber temp
    axes[0].plot(of_list, comb_temp, marker='o', linestyle='-')
    axes[0].set_ylabel('Chamber Temp (K)')
    axes[0].grid(True)
    axes[0].set_title('Chamber Temp (K) vs. O/F')

    # subplot 2, throat temp
    axes[1].plot(of_list, throat_temp, marker='o', linestyle='-')
    axes[1].set_ylabel('Throat Temp (K)')
    axes[1].grid(True)
    axes[1].set_title('Throat Temp (K) vs. O/F')

    # subplot 3, exit temp
    axes[2].plot(of_list, exhaust_temp, marker='o', linestyle='-')
    axes[2].set_ylabel('Exit Temp (K)')
    axes[2].grid(True)
    axes[2].set_title('Exit Temp (K) vs. O/F')

    # subplot 4, isp
    axes[3].plot(of_list, isp_list, marker='o', linestyle='-')
    axes[3].set_xlabel('O/F (MR)')
    axes[3].set_ylabel('Isp (s)')
    axes[3].grid(True)
    axes[3].set_title('Isp (s) vs. O/F')

    # show plot
    plt.tight_layout()
    plt.show()