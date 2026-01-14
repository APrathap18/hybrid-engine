def hoop_stress_calcs(r_c_m, t_m, pc_pa, y_s, pamb_pa=101325.0):
    """
    Finds the hoop stress on the walls of the chamber
    
    :param r_c_m: Radius of the chamber (m)
    :param t_m: Thickness of the chamber walls (m)
    :param pc_pa: Chamber pressure (Pa)
    :param y_s: Yield strength (Pa)
    """
    # Thin-walled cylinder assumption
    if t_m < r_c_m/10:
        hoop_stress = pc_pa * r_c_m / t_m # Hoop stress (Pa)
    # Thick-walled cylinder assumption (Lame equation)
    elif t_m >= r_c_m/10:
        r_o_m = r_c_m + t_m
        # Finds hoop stress at inner wall (Maximum)
        hoop_stress = (r_c_m**2 * pc_pa - r_o_m**2 * pamb_pa) / (r_o_m**2 - r_c_m**2) - (r_c_m**2 * r_o_m**2 * (pamb_pa - pc_pa)) / (r_c_m**2 * (r_o_m**2 - r_c_m**2))
    fos = y_s / hoop_stress

    print(f"Hoop Stress (Pa): {hoop_stress} Pa")
    print(f"Factor of Safety: {fos}")

    return hoop_stress