from rocketcea.cea_obj import add_new_fuel
import nozzle_sizing
import structural_calcs
import thermal_calcs

INCHES_TO_METERS = 0.0254
METERS_TO_INCHES = 39.3701

# ----------------------------
# Custom fuel: Paraffin/Al/C40H82, 80/10/10 by mass
# ----------------------------
paraffin_htpb_carbon_card = """
fuel Paraffin80     C 73.0    H 148.0               wt%=85.0
h,cal=-4.4464E+05   t(k)=298.15

fuel HTPB12         C 7.3165  H 10.3360  O 0.1063   wt%=12.0
h,cal=1200.0        t(k)=298.15          rho = 0.9220

fuel CarbonBlack3   C 1.0                           wt%=3.0
h,cal=0             t(k)=298.15
"""

add_new_fuel("Paraffin_HTPB_carbon_black", paraffin_htpb_carbon_card)

import matplotlib.pyplot as plt
import throat_sizing
import plot_OF
import orifice_sizing
import stress_calcs

# ----------------------------
# Globals / design parameters
# ----------------------------
oxName = 'N2O'
fuelName = 'Paraffin_HTPB_carbon_black'
pamb = 14.7 # psia
eps = 4

# Densities (design values)
rhoN2O  = 750.0  # kg/m^3, liquid N2O at ~20–25 C
rhofuel = 950.0  # kg/m^3, paraffin/Al/C40 mixture

def main():
    of = 6
    pc = 300 # psia
    F = 890 # N
    p1_pa = 360 * 6894.76 # Pa (300 psi chamber + 60 psi drop over injector)
    pc_pa = 300 * 6894.76 # Pa
    Cd = 0.7
    FS = 3
    burn_time = 5 # seconds

    # Chamber geometry (inches)
    L_star = 60.0   # in
    D_c    = 2.5   # in
    conv_angle = 30 # degrees
    div_angle = 15  # degrees

    # Material Properties
    therm_cond = 5
    density = 5
    c_p = 5
    
    [eps, Cf, A_t, A_e, Dt_m, De_m, m_dot, mach_exit, T0, cstar] = nozzle_sizing.nozzle_sizer(F, pc_pa, of, oxName, fuelName)

    At_in = A_t / (INCHES_TO_METERS**2)
    Dt_in = Dt_m / INCHES_TO_METERS
    De_in = De_m / INCHES_TO_METERS
    Dc_m = D_c * INCHES_TO_METERS

    # [At, At_in, m_dot, throat_dia, cstar, Tc_K] = throat_sizing.throat_sizing_function(of, pc, F, eps, pamb, rhoN2O, rhofuel, oxName, fuelName) # in^2
    plot_OF.plot_OF(pc, eps, oxName, fuelName, pamb)

    _, _, L_straight_in, L_conv_in, V_c_in3, L_div_in = orifice_sizing.orifice_area(At_in, of, Cd, m_dot, L_star, D_c, Dt_in, De_in, conv_angle, div_angle, p1_pa, pc_pa, rhoN2O)
    
    L_total = L_straight_in * INCHES_TO_METERS + L_conv_in * INCHES_TO_METERS + L_div_in * INCHES_TO_METERS

    fos = structural_calcs.hoop_stress_calcs(Dc_m/2, 0.018, pc_pa, 9.65e7)

    thermal_calcs.calculate_wall_temperature(300, 0.1, T0, 0.01, pc_pa, oxName, fuelName, of, eps, mach_exit, burn_time, Dc_m, cstar, Dt_m, De_m)
if __name__ == "__main__":
    main()