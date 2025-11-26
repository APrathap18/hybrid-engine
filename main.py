import numpy as np
from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import matplotlib.pyplot as plt
import throat_sizing
import plot_OF
import orifice_sizing

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
eps = 4

def main():
    of = 1.4
    pc = 300 # psia
    F = 890 # N

    At_in = throat_sizing.throat_sizing_function(of, pc, F) # in^2

    # EPS of 4 estimated
    plot_OF.plot_OF(pc, eps)

    orifice_sizing.orifice_area(At_in, of)
    #print(Dt_in)

if __name__ == "__main__":
    main()