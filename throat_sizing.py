import numpy as np
from rocketcea.cea_obj import CEA_Obj

m_sq_to_in_sq = 1550 #1 m^2 = 1550 in^2


# This is a function that takes O/F ratio, chamber pressure (psia), 
# and thrust (N) as inputs, and outputs the diameter (in) of the throat 
def throat_sizing_function(of, pc, F, eps, pamb, n2o_density, paraffin_density, oxname, fuelname):
    # creates CEA object with Ox and Fuel
    cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)

    # 1) Temperatures (RocketCEA temps are typically in degR; degR -> K = * 5/9)
    Tc_R, Tt_R, Te_R = cea_obj.get_Temperatures(Pc=pc, MR=of, eps=eps)
    Tc_K, Tt_K, Te_K = Tc_R * (5/9), Tt_R * (5/9), Te_R * (5/9)

    # 2) Use *ambient-corrected* thrust coefficient
    # Ambient-corrected accounts for ambient pressure pushing on the nozzle exit
    CF, CFamb, mode = cea_obj.get_PambCf(Pamb=pamb, Pc=pc, MR=of, eps=eps)

    # 3) Get cstar and convert from ft/s to m/s
    cstar = cea_obj.get_Cstar(Pc = pc, MR = of) * 0.3048

    # 4) Get ISP_ambient = (CFamb * C*) / g0
    isp = (cstar * CFamb) / 9.80665

    # get exhaust velocity
    v_e = isp * 9.80655 # m/s (assumed 90% of ideal) max - removed 0.9 , just including a
                        # general cstar efficiency

    # calculate mdot
    mdot = F/v_e #kg/s

    # calculate mdot for N2O and Paraffin
    n2o_mdot = (mdot * of) / (1 + of)
    fuel_mdot = mdot / (1 + of)

    # convert pc to pascals
    pc_pa = 6894.76 * pc #Pa
    # calculate area of the throat (m^2)
    At = mdot * cstar/pc_pa
    # calculate diameter of throat (m)
    Dt_m = 2 * np.sqrt(At/np.pi)
    # convert to inches
    Dt_in = Dt_m / 0.3048 * 12

    print("----------------------------------")
    print("Throat Sizing Parameters")
    print("----------------------------------")
    print(f"Throat Area (m^2): {At}")
    print(f"Throat Area (in^2): {At * m_sq_to_in_sq}")
    print(f"Throat Diameter (in): {Dt_in}")
    print("----------------------------------")
    print("Flow Parameters")
    print("----------------------------------")
    print(f"ISP (s): {isp}")
    print(f"C* (m/s): {cstar}")
    print(f"Exhaust Velocity (m/s): {v_e}")
    print(f"M Dot (kg/s): {mdot}")
    print(f"N2O M Dot (kg/s): {n2o_mdot}")
    print(f"Paraffin M Dot (kg/s): {fuel_mdot}")
    print(f"N2O Volumetric Flow Rate (m^3/s): {n2o_mdot/n2o_density}")
    print(f"N2O Volumetric Flow Rate (in^3/s): {61020 * n2o_mdot/n2o_density}")
    print(f"Paraffin Volumetric Flow Rate (m^3/s): {fuel_mdot/paraffin_density}")
    print(f"Paraffin Volumetric Flow Rate (in^3/s): {61020 * fuel_mdot/paraffin_density}")
    print("----------------------------------")
    print("Temperature Parameters")
    print("----------------------------------")
    print(f"Combustion Temperature (K): {Tc_K}")
    print(f"Throat Temperature (K): {Tt_K}")
    print(f"Exhaust Temperature (K): {Te_K}")
    print("----------------------------------")

    return(At, At * m_sq_to_in_sq, mdot, Dt_m, cstar, Tc_K)