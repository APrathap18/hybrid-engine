import numpy as np
from rocketcea.cea_obj import CEA_Obj
import matplotlib.pyplot as plt
import throat_sizing
import plot_OF
import TTP_Orifice
#import TTP Orifice

oxName = 'GOX'
fuelName = 'RP1'
pamb = 14.7 # psia

def main():
    of = 1.4
    pc = 150 # psia
    F = 15 # N

    At_in = throat_sizing.throat_sizing_function(of, pc, F) # in^2

    # EPS of 1 because no diverging section
    plot_OF.plot_OF(pc, 1)

    TTP_Orifice.orifice_area(At_in)
    #print(Dt_in)

if __name__ == "__main__":
    main()