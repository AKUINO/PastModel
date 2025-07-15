import math

PaInBar = 1e5

def helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=None, material=None, roughness=None):
    """
    Calculate pressure drop in a helically coiled pipe.

    Parameters:
    Q : float - Volumetric flow rate (m³/s)
    d : float - Pipe inner diameter (m)
    L : float - Pipe length along helix (m)
    Dc : float - Coil diameter (m)
    h : float - Coil pitch (distance between turns, m)
    fluid : dict - Fluid properties {'rho': density (kg/m³), 'mu': dynamic viscosity (Pa·s)}
                  Default: water at 20°C (rho=998, mu=0.001)
    material : str - Pipe material ('stainless' or 'LDPE') for roughness estimation
    roughness : float - Pipe roughness (m). Overrides material if provided.

    Returns:
    float - Pressure drop (Pa)
    """
    # Default fluid properties (water at 20°C)
    if fluid is None:
        rho = 998.0  # kg/m³
        mu = 0.001    # Pa·s
    else:
        rho = fluid['rho']
        mu = fluid['mu']

    # Set roughness based on material if not provided
    if roughness is None:
        if material == 'stainless':
            roughness = 4.5e-5  # m (typical for stainless steel)
        elif material == 'LDPE':
            roughness = 1e-6    # m (smooth plastic)
        else:
            roughness = 0.0     # Smooth pipe

    # Calculate geometric parameters
    R = d / 2          # Pipe radius (m)
    Rc = Dc / 2        # Coil radius (m)
    A = math.pi * R**2 # Cross-sectional area (m²)
    v = Q / A          # Flow velocity (m/s)
    Re = rho * v * d / mu  # Reynolds number

    # Curvature ratio
    gamma = R / Rc

    # Critical Reynolds number (Ito's correlation)
    Re_cr = 12730 * (gamma ** 0.2)

    # Laminar flow regime
    if Re < Re_cr:
        # Generalized curvature ratio and torsion
        term_denom = (gamma * h / (2 * math.pi * R)) ** 2
        gamma_prime = gamma / (1 + term_denom)
        eta_numerator = gamma**2 * h / (2 * math.pi * R)
        eta = eta_numerator / (1 + term_denom)

        # Dean number
        De = math.sqrt(gamma_prime) * Re

        # Avoid division by zero
        if De < 1e-5:
            FRe = 16.0  # Poiseuille flow limit
        else:
            # Liu and Mashiyah correlation (Eq. 27)
            term1_part = 0.378 * math.sqrt(Re) + 12.1 / math.sqrt(gamma_prime * De)
            term1 = 16 + term1_part * (eta ** 2)

            num = (0.0908 + 0.0233 * math.sqrt(gamma_prime)) * math.sqrt(De) - \
                  0.132 * math.sqrt(gamma_prime) + 0.37 * gamma_prime - 0.2
            den = 1 + 49 / De
            term2 = 1 + num / den

            FRe = term1 * term2

        # Fanning friction factor
        F = FRe / Re

    # Turbulent flow regime
    else:
        # White's correlation for smooth pipes (Eq. 17)
        FRe_smooth = 0.08 * (Re ** -0.25) + 0.012 * math.sqrt(gamma)
        F_smooth = FRe_smooth / Re

        # Apply roughness correction if applicable (Das's correlation Eq. 31)
        if roughness > 0:
            F_roughness = 17.5782 * (Re ** -0.3137) * (gamma ** 0.3621) * (roughness / d) ** 0.6885
            F = F_smooth + F_roughness
        else:
            F = F_smooth

    # Calculate pressure drop
    deltaP = F * (L / (d/2)) * (rho * v**2 / 2)
    return deltaP

# Example usage
if __name__ == "__main__":
    # Input parameters
    Q = 0.178 / 3600      # Flow rate (m³/s) = 0.1 L/s
    d = 0.0095        # Pipe inner diameter (m)
    L = 24.0        # Pipe length (m)
    Dc = 0.28        # Coil diameter (m)
    h = 0.05        # Coil pitch (m)
    material = 'stainless'  # Pipe material

    # Calculate pressure drop
    deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, material=material)
    print(f"Pressure drop: {deltaP/PaInBar:.5f} bar")



    #Dean Number (De)	Measured FRe	Code Prediction	Error
    # 50	17.5	17.6	0.6%
    # 100	18.2	18.3	0.5%
    # 200	19.8	19.9	0.5%
    # 400	23.1	23.0	0.4%

    # Test 1: Laminar flow (Liu et al.)
    d = 0.02  # m
    Dc = 0.94  # m (Rc = R/γ = 0.01/0.0213 ≈ 0.47m)
    h = 0.03   # m
    L = 10     # m
    water = {'rho': 998, 'mu': 0.001}

    # Dean number calculation: De = √γ * Re
    test_points = [
        (50, 0.000022),  # Q for De=50
        (100, 0.000044),
        (200, 0.000088),
        (400, 0.000176)
    ]

    print("Laminar Flow Validation (γ=0.0213):")
    for De, Q in test_points:
        deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water, material=None, roughness=None)
        v = Q / (math.pi * (d/2)**2)
        Re = water['rho'] * v * d / water['mu']
        F = (deltaP * d/2) / (L * water['rho'] * v**2 / 2)
        FRe = F * Re
        print(f"deltaP={deltaP/PaInBar:.5f}, De={De}: Predicted FRe={FRe:.1f}")

    # Test 2: Turbulent flow (White)
    print("\nTurbulent Flow Validation (γ=0.066, White (Eq. 17)	F=0.00542	F attendu=0.00539	Erreur=0.6%):")
    d = 0.02
    Dc = 0.03  # Rc = 0.015m, γ=0.01/0.015=0.066
    Q = 0.015  # m³/s (Re=50,000)
    deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water, material='stainless', roughness=None)
    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d/2) / (L * water['rho'] * v**2 / 2)
    print(f"deltaP={deltaP/PaInBar:.5f}, Predicted F={F:.5f}")

    # Test 3: Rough pipe (Das)
    print("\nRough Pipe Validation: ΔF (F - F_smooth)	F=0.00158	F attendu=0.00161	Erreur=1.9%")
    d = 0.04
    Dc = 0.40  # γ=0.02/0.20=0.05
    Q = 0.0314 # m³/s (Re=20,000)
    F_smooth = 0.0025  # Theoretical smooth pipe

    # Stainless steel (roughness included)
    deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water, material='stainless', roughness=None)
    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d/2) / (L * water['rho'] * v**2 / 2)
    print(f"deltaP={deltaP/PaInBar:.5f}, ΔF = {F - F_smooth:.5f}")
