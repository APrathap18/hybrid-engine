import numpy as np
from rocketcea.cea_obj import CEA_Obj
import warnings
import matplotlib.pyplot as plt
import traceback

def plot_OF(pc, eps, oxName, fuelName, pamb):
    #print("\n[DEBUG] plot_OF called from:")
    #traceback.print_stack(limit=4)

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
        with warnings.catch_warnings():
            warnings.simplefilter("error", RuntimeWarning)
            try:
                # 1) Pre-check gamma to make sure it isn't a bad gamma
                # Gets molar weight and gamma at exit
                mw_e, gam_e = cea_obj.get_exit_MolWt_gamma(Pc = pc, MR = of, eps = eps)

                # 2) Temperatures (RocketCEA temps are typically in degR; degR -> K = * 5/9)
                Tc_R, Tt_R, Te_R = cea_obj.get_Temperatures(Pc=pc, MR=of, eps=eps)
                Tc_K, Tt_K, Te_K = Tc_R * (5/9), Tt_R * (5/9), Te_R * (5/9)

                # 3) Use *ambient-corrected* thrust coefficient
                # Ambient-corrected accounts for ambient pressure pushing on the nozzle exit
                CF, CFamb, mode = cea_obj.get_PambCf(Pamb=pamb, Pc=pc, MR=of, eps=eps)

                # Makes sure the ambient CF is good
                if not np.isfinite(CFamb):
                    continue
                
                # 4) Get cstar and convert from ft/s to m/s
                cstar = cea_obj.get_Cstar(Pc = pc, MR = of) * 0.3048

                # 5) Get ISP_ambient = (CFamb * C*) / g0
                Isp = (cstar * CFamb) / 9.80665

                # Since everything is successful, append to list

                # add o/f ratio to list
                of_list.append(of)

                # append temps to corresponding lists, converting to Kelvin
                comb_temp.append(Tc_K)
                throat_temp.append(Tt_K)
                exhaust_temp.append(Te_K)

                # add isp to list
                isp_list.append(Isp)
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
    axes[0].set_ylabel('Chamber Temp [K]')
    axes[0].grid(True)
    axes[0].set_title('Chamber Temp [K] vs. O/F')

    # subplot 2, throat temp
    axes[1].plot(of_list, throat_temp, marker='o', linestyle='-')
    axes[1].set_ylabel('Throat Temp [K]')
    axes[1].grid(True)
    axes[1].set_title('Throat Temp [K] vs. O/F')

    # subplot 3, exit temp
    axes[2].plot(of_list, exhaust_temp, marker='o', linestyle='-')
    axes[2].set_ylabel('Exit Temp [K]')
    axes[2].grid(True)
    axes[2].set_title('Exit Temp [K] vs. O/F')

    # subplot 4, isp
    axes[3].plot(of_list, isp_list, marker='o', linestyle='-')
    axes[3].set_xlabel('O/F [MR]')
    axes[3].set_ylabel('Isp [s]')
    axes[3].grid(True)
    axes[3].set_title('Isp [s] vs. O/F')

    # show plot
    plt.tight_layout()
    plt.show()