from rocketcea.cea_obj import add_new_fuel

INCHES_TO_METERS = 0.0254

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

import matplotlib.pyplot as plt
import throat_sizing
import plot_OF
import orifice_sizing
import stress_calcs

# ----------------------------
# Globals / design parameters
# ----------------------------
oxName = 'N2O'
fuelName = 'Paraffin_Al_C40_80_10_10'
pamb = 14.7 # psia
eps = 4

# Densities (design values)
rhoN2O  = 750.0  # kg/m^3, liquid N2O at ~20–25 C
rhofuel = 950.0  # kg/m^3, paraffin/Al/C40 mixture

def main():
    of = 5
    pc = 300 # psia
    F = 890 # N
    p1 = 360 * 6894.76 # Pa (300 psi chamber + 60 psi drop over injector)
    p2 = 300 * 6894.76 # Pa
    Cd = 0.7
    FS = 3
    burn_time = 5 # seconds

    # Chamber geometry (inches)
    L_star = 45.0   # in
    D_c    = 1.25   # in
    L_c_in = 0.0    # in
    L_c_m = 0.0     # in
    V_c_in3 = 0.0   # in^3

    # Material Properties
    therm_cond = 5
    density = 5
    c_p = 5
    

    [At_in, m_dot, throat_dia_m, cstar, Tc_K] = throat_sizing.throat_sizing_function(of, pc, F, eps, pamb, rhoN2O, rhofuel, oxName, fuelName) # in^2

    # EPS of 4 estimated
    plot_OF.plot_OF(pc, eps, oxName, fuelName, pamb)

    _, _, L_c_in, V_c_in3 = orifice_sizing.orifice_area(At_in, of, Cd, m_dot, L_star, D_c, p1, p2, rhoN2O, rhofuel, pamb, oxName, fuelName)
    #print(Dt_in)

    #L_c_m = L_c_in * 0.0254

    #thick = []
    #h_g = []

    #thick, h_g, _, _ = stress_calcs.wall_thickness(pc, of, eps, cstar, throat_dia_m, D_c * INCHES_TO_METERS, Tc_K, burn_time, therm_cond, density, c_p, T_table, yield_strength, oxName, fuelName, FS, L_c_m, p_amb=101325.0, t_max0=0.05)
if __name__ == "__main__":
    main()