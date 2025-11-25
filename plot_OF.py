import numpy as np
from rocketcea.cea_obj import CEA_Obj
import matplotlib.pyplot as plt
oxName = 'GOX'
fuelName = 'RP1'
pamb = 14.7 # psia
# This is a test

def plot_OF(pc, eps):
    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxName, fuelName = fuelName)

    # lists for chamber temp, throat temp, exhaust temp, isp, and o/f ratios
    comb_temp = []
    throat_temp = []
    exhaust_temp = []
    isp_list = []

    of_list = []

    # iterates through o/f ratios in 0.5 step sizes
    for of in np.arange(0.5, 10.0, 0.25):
        try:
            # add o/f ratio to list
            of_list.append(of)

            # get all temperatures in one list
            temp = cea_obj.get_Temperatures(Pc = pc, MR = of, eps = eps)

            # append temps to corresponding lists, converting to Kelvin
            comb_temp.append(temp[0] * 5/9)
            throat_temp.append(temp[1] * 5/9)
            exhaust_temp.append(temp[2] * 5/9)

            # get fuel coefficient
            Cf = cea_obj.get_PambCf(Pamb = 14.7, Pc = pc, MR = of, eps = eps)[0]

            # get cstar
            cstar = cea_obj.get_Cstar(Pc = pc, MR = of) * 0.3048

            # add isp to list
            isp_list.append(cstar * Cf / 9.80655)
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