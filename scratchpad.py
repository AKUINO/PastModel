def calcul_performance():
    # Vitesses
    v_tube = Q / (np.pi * d_int**2 / 4)
    v_ann = Q / (np.pi * (d_calandre**2 - d_ext**2) / 4)

    # Reynolds (crème à 70°C)
    ρ = 1020  # kg/m³
    μ = 0.001  # Pa.s
    Re_tube = ρ * v_tube * d_int / μ
    Re_ann = ρ * v_ann * (d_calandre - d_ext) / μ

    # Coefficients convection (Gnielinski)
    Nu_tube = 0.023 * Re_tube**0.8 * 4.5**0.4
    h_tube = Nu_tube * 0.5 / d_int

    # [...] (calculs similaires pour annulaire)

    # U global
    U = 1 / (1/h_tube + 1/h_ann + 1e-4)

    # Puissance
    ΔT = 21.3  # K
    P = U * (np.pi * d_ext * L) * ΔT

    return P, Re_tube, Re_ann