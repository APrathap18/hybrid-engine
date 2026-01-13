def hoop_stress_calcs(r_c, t, pc, y_s):
    hoop_stress = pc * r_c / t # Hoop stress (Pa)
    fos = y_s / hoop_stress

    return [hoop_stress, fos]