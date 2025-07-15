import math

def perte_charge_simple(Re, D=0.010, Rc=0.15, L=10.0):
    rho = 998
    mu = 1.002e-3
    nu = mu / rho

    v = Re * nu / D
    A = math.pi * D**2 / 4
    Q = v * A * 3600  # m³/h

    gamma = D / (2 * Rc)
    De = Re * math.sqrt(gamma)

    if De < 200:
        f_Fanning = 0.47136 * De**0.25
    else:
        f_Fanning = 0.556 + 0.0969 * math.sqrt(De)

    f_Darcy = 4 * f_Fanning

    print(f"Re = {Re}, Q = {Q:.4f} m³/h, v = {v:.4f} m/s, De = {De:.1f}, f_Darcy = {f_Darcy:.4f}")
    return f_Darcy

# Boucle sur les Re testés
for Re in [500, 1000, 1500, 2000, 3000]:
    perte_charge_simple(Re)