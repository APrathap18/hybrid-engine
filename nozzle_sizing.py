import math
from rocketcea.cea_obj import CEA_Obj

PSIA_TO_PA = 6894.757293168
PA_TO_PSIA = 1.0 / PSIA_TO_PA
FT_TO_M = 0.3048
m_sq_to_in_sq = 1550 #1 m^2 = 1550 in^2

def _scalar(x):
    """
    Checks if x is a scalar. If it isn't (i.e. tuple or list), it returns the first element
    
    :param x: Data of some type 
    """
    return x[0] if isinstance(x, (tuple, list)) else x

def pc_over_pe_from_eps(cea, pc_psia, of, eps):
    """
    Solves for a value epsilon = pc/pe
    The ideal is that exhaust pressure = ambient pressure (perfectly expanded)
    The best case is when eps is passed to this function so that pc/pe = pc/pamb
    
    :param cea: cea object
    :param pc_psia: chamber pressure (psia)
    :param of: of ratio
    :param eps: expansion ratio
    """
    pc_over_pe = cea.get_PcOvPe(Pc = pc_psia, MR = of, eps = eps)
    return _scalar(pc_over_pe)

def find_eps_for_pe_equals_pamb(cea, pc_psia, of, pamb_psia, eps_low = 1.01, eps_high = 200.0, iters = 70):
    """
    Finds the expansion ratio for which the exit pressure equals the ambient pressure

    :param cea: NASA CEA object
    :param pc_psia: chamber pressure (psia)
    :param of: of ratio
    :param pamb_psia: ambient pressure (psia)
    :param eps_low: lower bound for eps
    :param eps_high: higher bound for eps
    :param iters: number of iterations
    """
    # Ideal value for pc/pe (pe = pamb)
    target = pc_psia / pamb_psia

    # Low bound for the value of pc/pe (when eps is very low)
    # This will give us a negative value for pc since target > the result from CEA
    # That's good, we want to find the point where f = 0 because then that means target = the pc/pe from CEA
    f_low = pc_over_pe_from_eps(cea, pc_psia, of, eps_low) - target

    # High bound for the value of pc/pe (when eps is very high)
    # This value will be positive
    f_high = pc_over_pe_from_eps(cea, pc_psia, of, eps_high) - target

    # This is essentially a binary search
    for _ in range(iters):
        # Average of low and high eps
        eps_mid = 0.5 * (eps_low + eps_high)
        
        # The value of pc/pe at this eps
        f_mid = pc_over_pe_from_eps(cea, pc_psia, of, eps_mid) - target

        # If the product is negative, that means that either f_low or f_mid is negative
        # That means that this section of values contains the root (where f(epsilon) = 0)
        # Therefore, change the max value to mid
        if f_low * f_mid <= 0:
            f_high = f_mid
            eps_high = eps_mid
        # If the product is positive, that means that either both are negative or both are positive
        # That means that the sign change is between mid and high, so update low accordingly
        else:
            eps_low = eps_mid
            f_low = f_mid

    return (0.5 * (eps_low + eps_high))

    
def nozzle_sizer(F_N, pc_pa, of, oxname, fuelname, pamb_pa=101325.0):
    """
    Puts all helper functions together with sizing functions to find the size of each nozzle section

    :param F_N: thrust (N)
    :param pc_pa: chamber pressure (pa)
    :param of: of ratio
    :param oxname: name of oxidizer (string)
    :param fuelname: name of fuel (string)
    :param pamb_pa: ambient pressure (pa)
    """
    cea = CEA_Obj(oxName = oxname, fuelName = fuelname)

    pc_psia = pc_pa * PA_TO_PSIA
    pamb_psia = pamb_pa * PA_TO_PSIA

    # 1. Find eps such that pe = pamb
    eps = find_eps_for_pe_equals_pamb(cea, pc_psia, of, pamb_psia)

    # 2. Find the fuel coefficient at ambient pressure for that eps
    Cf = cea.get_PambCf(Pamb = pamb_psia, Pc = pc_psia, MR = of, eps = eps)
    Cf = _scalar(Cf)

    # 3. Find throat area from the equation F = Cf * Pc * At
    A_t = F_N / (Cf * pc_pa)

    # 4. Find exit area using the eps
    A_e = A_t * eps

    # Find diameters of throat and exit
    D_t = math.sqrt(4.0 * A_t / math.pi)
    D_e = math.sqrt(4.0 * A_e / math.pi)

    # Temperatures (RocketCEA temps are typically in degR; degR -> K = * 5/9)
    Tc_R, Tt_R, Te_R = cea.get_Temperatures(Pc=pc_psia, MR=of, eps=eps)
    Tc_K, Tt_K, Te_K = Tc_R * (5/9), Tt_R * (5/9), Te_R * (5/9)

    # Get cstar and convert from ft/s to m/s
    cstar = cea.get_Cstar(Pc = pc_psia, MR = of) * 0.3048
    
    # Get ISP_ambient = (CFamb * C*) / g0
    isp = (cstar * Cf) / 9.80665

    # Find molar weight and gamma
    mw, gamma = cea.get_Chamber_MolWt_gamma(Pc=pc_psia, MR=of, eps=eps)
    
    # Gas constant
    RU = 8314.462618  # J/(kmol*K)
    Rspec = RU / mw   # J/(kg*K)

    # Local area of sound at throat (Mach is 1 so this is also the velocity at the throat)
    a_t = math.sqrt(gamma * Rspec * Tt_K)

    # get exhaust velocity
    v_e = isp * 9.80655

    # Get mass flow rate (kg/s)
    m_dot = pc_pa * A_t / cstar 

    print("----------------------------------")
    print("Nozzle Sizing Parameters")
    print("----------------------------------")
    print(f"Throat Area (m^2): {A_t} m^2")
    print(f"Throat Area (in^2): {A_t * m_sq_to_in_sq} in^2")
    print(f"Throat Diameter (in): {D_t / 0.3048 * 12} in")
    print(f"Nozzle Exit Area (m^2): {A_e} m^2")
    print(f"Nozzle Exit Area (in^2): {A_e * m_sq_to_in_sq} in^2")
    print(f"Nozzle Exit Diameter (in): {D_e / 0.3048 * 12} in")
    print(f"Nozzle Expansion Ratio: {eps}")
    print("----------------------------------")
    print("Flow Parameters")
    print("----------------------------------")
    print(f"Propellant Mass Flow Rate (kg/s): {m_dot} kg/s")
    print(f"Chamber Temperature (K): {Tc_K} K")
    print(f"Throat Temperature (K): {Tt_K} K")
    print(f"Exhaust Temperature (K): {Te_K} K")
    print(f"Throat Velocity (m/s): {a_t} m/s")
    print(f"Exit Velocity (m/s): {v_e} m/s")

    return(eps, Cf, A_t, A_e, D_t, D_e, m_dot)

# def nozzle_sizing(F, pc, pe, of, gamma, oxname, fuelname):
#     # creates CEA object with Ox and Fuel
#     cea_obj = CEA_Obj(oxName = oxname, fuelName = fuelname)

#     pc_psia = pc * PA_TO_PSIA
#     pe_psia = pe * PA_TO_PSIA

#     # Temperatures (RocketCEA temps are typically in degR; degR -> K = * 5/9)
#     T0_R, _, _ = cea_obj.get_Temperatures(Pc=pc, MR=of, eps=1)
#     T0 = T0_R * (5/9)

#     # 1. Critical (Throat) Temperature (K)
#     T_t = 2 * T0 / (gamma + 1)

#     # 2. Throat velocity (local speed of sound)
#     a = math.sqrt(gamma * R * T_t)

#     # 3. Exit velocity (m/s)
#     V = math.sqrt((2 * gamma) / (gamma - 1) * R * T0 * (1 - (pe / pc)**((gamma-1)/gamma)))

#     # 4. Propellant Mass Flow Rate (kg/s)
#     m_dot = F / a

#     # 5. Specific Volume at Nozzle Entrance (m^3/kg)
#     v_1 = R * T0 / pc

#     # 6. Throat Specific Volume (m^3/kg)
#     v_t = v_1 * ((gamma + 1) / 2)**(1 / (gamma-1))

#     # 7. Nozzle Exit Specific Volume (m^3/kg)
#     v_e = (10**((gamma-1)/gamma) * (pc/pe)) / (10**(gamma-1) * v_1)

#     # 8. Throat Area (m^2)
#     A_t = m_dot * v_t / a

#     # 9. Nozzle Exit Area (m^2)
#     A_e = m_dot * v_e / V
    
#     # 10. Throat Diameter (m)
#     Dt_m = 2 * math.sqrt(A_t/math.pi)

#     # 11. Exit Diameter (m)
#     Dt_e = 2 * math.sqrt(A_e/math.pi)

#     print("----------------------------------")
#     print("Nozzle Sizing Parameters")
#     print("----------------------------------")
#     print(f"Throat Area (m^2): {A_t} m^2")
#     print(f"Throat Area (in^2): {A_t * m_sq_to_in_sq} in^2")
#     print(f"Throat Diameter (in): {Dt_m / 0.3048 * 12} in")
#     print(f"Nozzle Exit Area (m^2): {A_e} m^2")
#     print(f"Nozzle Exit Area (in^2): {A_e * m_sq_to_in_sq} in^2")
#     print(f"Nozzle Exit Diameter (in): {Dt_e / 0.3048 * 12} in")
#     print(f"Nozzle Expansion Ratio: {A_e/A_t}")
#     print("----------------------------------")
#     print("Flow Parameters")
#     print("----------------------------------")
#     print(f"Propellant Mass Flow Rate (kg/s): {m_dot} kg/s")
#     print(f"Throat Temperature (K): {T_t} K")
#     print(f"Throat Velocity (m/s): {a} m/s")
#     print(f"Exit Velocity (m/s): {V} m/s")
#     print(f"Specific Volume at Nozzle Entrance (m^3/kg): {v_1} m^3/kg")
#     print(f"Specific Volume at Throat (m^3/kg): {v_t} m^3/kg")
#     print(f"Specific Volume at Nozzle Exit (m^3/kg): {v_e} m^3/kg")

