# Mesure à Q = 315 L/h → ΔP_mesuré = 1.4 bar (après correction)
Q = 315 / 3600  # Conversion L/h → m³/s
rho = 1000       # kg/m³ (eau)
L = 30           # m (longueur totale estimée)
D = 0.011        # m (diamètre intérieur)

k_perte = (2.5 * D**5) / (rho * Q**2 * L)
print (k_perte)

import numpy as np

def pertes_singulieres(Q, obstacles, diametres, rho=1000):
    """
    Calcule les pertes singulières dynamiques
    Args:
        Q: Débit (m³/s)
        obstacles: Liste de tuples [('type', geo_params), ...]
                   ex: [('coudes', 0.4, D_ref), ('retrecissement', D1, D2)]
        diametres: Dict des diamètres {nom_section: D(m)}
    Returns:
        ΔP_sing (bar)
    """
    ΔP = 0
    for obs in obstacles:
        type_obs = obs[0]

        if type_obs == 'coudes':
            K, D_ref, n = obs[1], obs[2], obs[3]
            A = np.pi * (D_ref/2)**2
            v = Q / A
            ΔP += n * K * (rho * v**2) / 2

        elif type_obs == 'retrecissement':
            D1, D2, n = obs[1], obs[2], obs[3]
            beta = min(D1, D2) / max(D1, D2)
            K = 0.5 * (1 - beta**2)
            A = np.pi * (min(D1, D2)/2)**2  # Section étroite
            v = Q / A
            ΔP += n * K * (rho * v**2) / 2

        elif type_obs == 'elargissement':
            D1, D2, n = obs[1], obs[2], obs[3]
            beta = min(D1, D2) / max(D1, D2)
            K = (1 - beta**2)**2
            A = np.pi * (min(D1, D2)/2)**2  # Section étroite
            v = Q / A
            ΔP += n * K * (rho * v**2) / 2

    return ΔP / 1e5  # Conversion Pa → bar

# Définir les obstacles
obstacles = [
    ('coudes', 0.4, 0.014, 6),      # 6 coudes dans section Ø14mm
    ('retrecissement', 0.0095, 0.0087, 18),  # 18 rétrécissements 9.5→8.7mm
    ('elargissement', 0.0087, 0.0095, 18)    # 18 élargissements 8.7→9.5mm
]

# Paramètres
Q = 145 / 3600 / 1000  # 145 L/h → m³/s
diametres = {'calandre': 0.014, 'jonction': 0.0095}

# Calcul
ΔP_sing = pertes_singulieres(Q, obstacles, diametres)
print(f"Pertes singulières totales: {ΔP_sing:.4f} bar")

def perte_charge_raccords(debit_l_h, d_tuyau_mm, d_retrec_mm, n_retrec, n_coudes):
    """
    Calcule la perte de charge totale due aux raccords (rétrécissements et coudes)

    Args:
        debit_l_h (float): Débit en litres par heure
        d_tuyau_mm (float): Diamètre principal du tuyau en mm
        d_retrec_mm (float): Diamètre au rétrécissement en mm
        n_retrec (int): Nombre de raccords avec rétrécissement
        n_coudes (int): Nombre de coudes à 90°

    Returns:
        tuple: (pertes_retrec_mCE, pertes_coudes_mCE, pertes_totales_mCE)
    """
    # Conversion des unités
    debit_m3_s = debit_l_h / 1000 / 3600  # Conversion L/h -> m³/s
    d_tuyau_m = d_tuyau_mm / 1000         # mm -> m
    d_retrec_m = d_retrec_mm / 1000       # mm -> m

    # Calcul des sections transversales
    section_tuyau = 3.1416 * (d_tuyau_m/2)**2
    section_retrec = 3.1416 * (d_retrec_m/2)**2

    # Calcul des vitesses
    v_tuyau = debit_m3_s / section_tuyau
    v_retrec = debit_m3_s / section_retrec

    # Coefficients de perte de charge (K)
    beta = d_retrec_m / d_tuyau_m
    K_retrec = 0.5 * (1 - beta**2)  # Coefficient pour rétrécissement brusque

    # Valeur validée pour les coudes à 90° standards
    K_coude = 0.5

    # Calcul des pertes (formule ΔP = K * (ρ * v²)/2)
    rho_eau = 1000  # kg/m³

    # Pertes pour les rétrécissements (basées sur v_retrec)
    pertes_retrec_pa = n_retrec * K_retrec * (rho_eau * v_retrec**2) / 2

    # Pertes pour les coudes (basées sur v_tuyau)
    pertes_coudes_pa = n_coudes * K_coude * (rho_eau * v_tuyau**2) / 2

    # Conversion Pa -> mCE (metres de colonne d'eau)
    pertes_retrec_mCE = pertes_retrec_pa / 100000
    pertes_coudes_mCE = pertes_coudes_pa / 100000
    pertes_totales_mCE = pertes_retrec_mCE + pertes_coudes_mCE

    return (pertes_retrec_mCE, pertes_coudes_mCE, pertes_totales_mCE)

# Application à votre cas
resultats_145 = perte_charge_raccords(
    debit_l_h=145,
    d_tuyau_mm=9.5,
    d_retrec_mm=8.7,
    n_retrec=18,
    n_coudes=6
)

resultats_310 = perte_charge_raccords(
    debit_l_h=310,
    d_tuyau_mm=9.5,
    d_retrec_mm=8.7,
    n_retrec=18,
    n_coudes=6
)

print(f"Pour 145 L/h : Pertes totales = {resultats_145[2]:.6f} bar")
print(f"Pour 310 L/h : Pertes totales = {resultats_310[2]:.6f} bar")