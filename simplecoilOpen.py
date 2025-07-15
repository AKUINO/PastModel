import math
# Affichage

PaInBar = 1e5


def perte_charge_serpentin(
        materiau: str,
        longueur_m: float,
        diametre_tube_m: float,
        debit_m3h: float,
        diametre_serpentin_m: float,
        pas_vertical_spires_m: float = None,
        temperature_C: float = 20.0
) -> dict:
    """
    Calcule la perte de charge d’un écoulement d’eau dans un serpentin hélicoïdal
    selon des corrélations issues de la littérature (Van Dyke, Hasson).

    Paramètres :
    - materiau : 'inox' ou 'lldpe' (à usage ultérieur pour la rugosité)
    - longueur_m : longueur développée du tube (m)
    - diametre_tube_m : diamètre intérieur du tube (m)
    - debit_m3h : débit volumique d'eau en m³/h
    - diametre_serpentin_m : diamètre du serpentin (m)
    - pas_vertical_spires_m : optionnel, actuellement ignoré
    - temperature_C : température de l’eau (°C), valeur par défaut 20°C

    Retourne :
    - Re : nombre de Reynolds
    - De : nombre de Dean
    - f_Fanning : facteur de friction Fanning
    - f_Darcy : facteur de friction Darcy-Weisbach
    - vitesse_m_s : vitesse moyenne de l'eau
    - perte_charge_Pa : perte de charge en Pascals
    """

    # Propriétés de l’eau à 20°C
    rho = 998       # kg/m³
    mu = 1.002e-3   # Pa.s
    nu = mu / rho   # m²/s

    # Calculs
    Q_m3s = debit_m3h / 3600  # conversion en m³/s
    A = math.pi * (diametre_tube_m**2) / 4
    v = Q_m3s / A
    Re = v * diametre_tube_m / nu
    gamma = diametre_tube_m / (2 * (diametre_serpentin_m / 2))  # D / (2 * Rc)
    De = Re * math.sqrt(gamma)

    # Facteur de friction (Van Dyke pour De<200, sinon Hasson)
    if De < 200:
        f_Fanning = 0.47136 * De**0.25
    else:
        f_Fanning = 0.556 + 0.0969 * math.sqrt(De)

    f_Darcy = 4 * f_Fanning
    delta_P = f_Darcy * (longueur_m / diametre_tube_m) * (rho * v**2 / 2)

    return {
        "Re": Re,
        "De": De,
        "f_Fanning": f_Fanning,
        "f_Darcy": f_Darcy,
        "vitesse_m_s": v,
        "perte_charge_Pa": delta_P
    }


def validate_code():
    # Fluid properties (water at 20°C)
    water = {'rho': 998, 'mu': 0.001}
    d = 0.02  # Pipe diameter (m)
    L = 10    # Pipe length (m)

    # Test 1: Laminar flow (Liu et al. 1994 data)
    print("Laminar Flow Validation (γ=0.0213):")
    Dc = 0.94  # Coil diameter (m) γ = R/Rc = 0.01/0.47 = 0.0213
    h = 0.03   # Coil pitch (m)

    # Corrected flow rates for target Dean numbers
    gamma = 0.01/0.47  # 0.0212766
    De_values = [50, 100, 200, 400]
    Q_values = [De / math.sqrt(gamma) * math.pi * d * water['mu'] / (4 * water['rho']) for De in De_values]

    for De, Q in zip(De_values, Q_values):
        deltaP = perte_charge_serpentin("inox", L, d, Q*3600, Dc, 0.05, 20)["perte_charge_Pa"]

        # Recalculate FRe from output
        A = math.pi * (d/2)**2
        v = Q / A
        Re = water['rho'] * v * d / water['mu']
        F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
        FRe = F * Re
        print(f"De={De}: deltaP={deltaP/100000:.5f} bar, Predicted FRe={FRe:.1f}")

    # Test 2: Turbulent flow (White 1929 data)
    print("\nTurbulent Flow Validation (γ=0.066):")
    d = 0.02
    Dc = 0.03  # γ = 0.01/0.015 = 0.0667
    h = 0.03
    Re_target = 50000
    Q = Re_target * math.pi * d * water['mu'] / (4 * water['rho'])  # Correct flow calculation
    deltaP = perte_charge_serpentin("inox", L, d, Q*3600, Dc, 0.05, 20)["perte_charge_Pa"]

    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
    print(f"deltaP={deltaP/100000:.5f} bar, Predicted F={F:.5f}")

    # Test 3: Rough pipe (Das 1993 data)
    print("\nRough Pipe Validation:")
    d = 0.04
    Dc = 0.80  # γ=0.02/0.40=0.05
    h = 0.03
    Re_target = 20000
    Q = Re_target * math.pi * d * water['mu'] / (4 * water['rho'])
    deltaP = perte_charge_serpentin("inox", L, d, Q*3600, Dc, 0.05, 20)["perte_charge_Pa"]

    v = Q / (math.pi * (d/2)**2)
    F = (deltaP * d) / (4 * L * water['rho'] * v**2 / 2)
    F_smooth = 0.0025  # Theoretical smooth pipe
    print(f"deltaP={deltaP/100000:.5f} bar, ΔF = {F - F_smooth:.5f}")

def valider_contre_experience():

    # Données expérimentales de De Amicis et al. (Figure 6)
    data_exp = [
        {"Re": 500, "f_Darcy_mesure": 0.31},
        {"Re": 1000, "f_Darcy_mesure": 0.22},
        {"Re": 1500, "f_Darcy_mesure": 0.17},
        {"Re": 2000, "f_Darcy_mesure": 0.145},
        {"Re": 3000, "f_Darcy_mesure": 0.12},
    ]

    results = []

    # Paramètres géométriques constants du test
    D = 0.010  # m (10 mm)
    Rc = 0.15  # m (300 mm de diamètre serpentin)
    L = 10.0   # m
    A = math.pi * D**2 / 4

    for entry in data_exp:
        Re = entry["Re"]
        # v = Re * nu / D
        nu = 1.002e-3 / 998  # viscosité cinématique
        v = Re * nu / D
        Q = v * A * 3600  # m³/h

        res = perte_charge_serpentin(
            materiau="inox",
            longueur_m=L,
            diametre_tube_m=D,
            debit_m3h=Q,
            diametre_serpentin_m=2*Rc
        )
        results.append({
            "Re": Re,
            "f_Darcy_calcule": res["f_Darcy"],
            "f_Darcy_mesure": entry["f_Darcy_mesure"]
        })
        print (f"{Re:.1f}\t{res['f_Darcy']:.5f}")
              # f"\t{res['perte_charge_Pa']/PaInBar:.3f}")



# Run validation
#validate_code()
valider_contre_experience()
