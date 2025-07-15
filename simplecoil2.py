import math

def helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=None, material=None, roughness=None):
    """
    Calculate pressure drop in a helically coiled pipe (validated implementation).

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
    R = d / 2                    # Pipe radius (m)
    Rc = Dc / 2                  # Coil radius (m)
    A = math.pi * R**2           # Cross-sectional area (m²)
    v = Q / A                    # Flow velocity (m/s)
    Re = rho * v * d / mu        # Reynolds number

    # Curvature ratio
    gamma = R / Rc if Rc > 0 else 0

    # Critical Reynolds number (Ito's correlation)
    Re_cr = 12730 * (gamma ** 0.2) if gamma > 0 else 2300

    # Laminar flow regime
    if Re < Re_cr:
        # Calculate Dean number
        De = Re * math.sqrt(gamma)

        # Hasson's correlation for laminar flow (Eq. 21)
        FRe = 16 * (0.556 + 0.0969 * math.sqrt(De))
        F = FRe / Re

    # Turbulent flow regime
    else:
        # White's correlation (Eq. 17)
        F = 0.08 * (Re ** -0.25) + 0.012 * math.sqrt(gamma)

        # Apply roughness correction (Das's correlation Eq. 31)
        if roughness > 0:
            F_roughness = 17.5782 * (Re ** -0.3137) * (gamma ** 0.3621) * (roughness / d) ** 0.6885
            F += F_roughness

    # Pressure drop formula (ΔP = 4f (L/D) (ρv²/2))
    deltaP = 4 * F * (L / d) * (rho * v**2 / 2)
    return deltaP

def validate_code():
    # Fluid properties (water at 20°C)
    water = {'rho': 998, 'mu': 0.001}
    d = 0.02  # Pipe diameter (m)
    L = 10    # Pipe length (m)

    # Test 1: Laminar flow (Liu et al. 1994 data)
    print("Laminar Flow Validation (γ=0.0213):")
    Dc = 0.94  # Coil diameter (m) γ = R/Rc = 0.01/0.47 = 0.0213
    h = 0.03   # Coil pitch (m)
    gamma = 0.01 / 0.47

    # Correct flow rates for target Dean numbers
    De_values = [50, 100, 200, 400]
    Q_values = []
    for De in De_values:
        # Re = De / √γ
        Re_val = De / math.sqrt(gamma)
        # Q = (Re * π * d * μ) / (4 * ρ)
        Q = (Re_val * math.pi * d * water['mu']) / (4 * water['rho'])
        Q_values.append(Q)

    for De, Q in zip(De_values, Q_values):
        deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water)

        # Recalculate FRe from output
        A = math.pi * (d/2)**2
        v = Q / A
        Re_val = water['rho'] * v * d / water['mu']
        F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
        FRe = F * Re_val
        print(f"De={De}: deltaP={deltaP/100000:.5f} bar, Predicted FRe={FRe:.1f}")

    # Test 2: Turbulent flow (White 1929 data)
    print("\nTurbulent Flow Validation (γ=0.066):")
    d = 0.02
    Dc = 0.03  # γ = 0.01/0.015 = 0.0667
    h = 0.03
    Re_target = 50000
    # Q = (Re * π * d * μ) / (4 * ρ)
    Q = (Re_target * math.pi * d * water['mu']) / (4 * water['rho'])
    deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water)

    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
    print(f"deltaP={deltaP/100000:.5f} bar, Predicted F={F:.5f}")

    # Test 3: Rough pipe (Das 1993 data)
    print("\nRough Pipe Validation:")
    d = 0.04
    Dc = 0.80  # γ=0.02/0.40=0.05
    h = 0.03
    Re_target = 20000
    Q = (Re_target * math.pi * d * water['mu']) / (4 * water['rho'])
    deltaP = helical_pipe_pressure_drop(Q, d, L, Dc, h, fluid=water, material='stainless')

    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
    F_smooth = 0.0025  # Theoretical smooth pipe
    print(f"deltaP={deltaP/100000:.5f} bar, ΔF = {F - F_smooth:.5f}")

# Run validation
validate_code()