import math
import numpy as np

# =============================================
# FONCTIONS POUR LES PROPRIÉTÉS PHYSIQUES DE L'EAU
# (Source : "Fundamentals of Heat and Mass Transfer" - Incropera & DeWitt)
# =============================================

def water_properties(T):
    """Retourne les propriétés de l'eau à température T (°C)"""
    T_kelvin = T + 273.15

    # Valeurs interpolées pour 20°C < T < 100°C
    rho = 1000 * (1 - (T_kelvin - 277.15)**2 / 508929.2)  # [kg/m³]
    k = 0.561 + 0.0029 * T - 1.01e-5 * T**2  # Conductivité [W/m·K]
    cp = 4186 - 1.15 * T + 0.004 * T**2  # Capacité thermique [J/kg·K]
    mu = 0.001 * 10**((1301/(998.333+8.1855*(T-20)+0.00585*(T-20)**2) - 3.30233))  # Viscosité [Pa·s]
    nu = mu / rho  # Viscosité cinématique [m²/s]
    alpha = k / (rho * cp)  # Diffusivité thermique [m²/s]
    Pr = nu / alpha  # Nombre de Prandtl [-]
    beta = 0.0002 * np.exp(0.005 * T)  # Coefficient de dilatation [K⁻¹]

    return {
        'rho': rho,    # Masse volumique [kg/m³]
        'k': k,        # Conductivité thermique [W/m·K]
        'cp': cp,      # Capacité calorifique [J/kg·K]
        'mu': mu,      # Viscosité dynamique [Pa·s]
        'nu': nu,      # Viscosité cinématique [m²/s]
        'Pr': Pr,      # Nombre de Prandtl [-]
        'beta': beta   # Coefficient de dilatation [K⁻¹]
    }

# =============================================
# CORRÉLATIONS DE TRANSFERT THERMIQUE
# =============================================

def h_internal(m_dot, Di, T_water, L=None):
    """Coefficient de convection interne (écoulement dans le tube)"""
    props = water_properties(T_water)
    Re = 4 * m_dot / (math.pi * Di * props['mu'])

    # Sélection de la corrélation selon le régime d'écoulement
    if Re < 2300:  # Laminaire
        if L is None:
            raise ValueError("Longueur du tube requise pour régime laminaire")
        Nu = 3.66 + (0.0668 * (Di/L) * Re * props['Pr']) / (1 + 0.04 * ((Di/L) * Re * props['Pr'])**(2/3))
    else:  # Turbulent (Dittus-Boelter)
        Nu = 0.023 * Re**0.8 * props['Pr']**0.4

    return Nu * props['k'] / Di

def h_external(Do, T_wall, T_bulk):
    """Coefficient de convection naturelle externe (Churchill-Chu)"""
    T_mean = 0.5 * (T_wall + T_bulk)
    props = water_properties(T_mean)
    delta_T = abs(T_wall - T_bulk)

    # Éviter delta_T = 0 pour les calculs
    if delta_T < 0.1:
        delta_T = 0.1

    Gr = (9.81 * props['beta'] * delta_T * Do**3) / (props['nu']**2)
    Ra = Gr * props['Pr']

    # Corrélation de Churchill-Chu
    term = (1 + (0.559 / props['Pr'])**(9/16))**(8/27)
    Nu = (0.60 + 0.387 * Ra**(1/6) / term)**2
    return Nu * props['k'] / Do

# =============================================
# CALCUL DU COEFFICIENT GLOBAL U
# =============================================

def calculate_U(params, tol=0.1, max_iter=50):
    """Calcule U par itération sur la température de paroi"""
    # Extraction des paramètres
    Di = params['Di']  # Diamètre intérieur [m]
    Do = params['Do']  # Diamètre extérieur [m]
    k_wall = params['k_wall']  # Conductivité inox [W/m·K]
    m_dot = params['m_dot']  # Débit massique [kg/s]
    T_coil_avg = params['T_coil_avg']  # Temp. moyenne eau serpentin [°C]
    T_tank = params['T_tank']  # Temp. eau cuve [°C]
    L = params.get('L', 10.0)  # Longueur du tube (valeur par défaut 10m)

    # Rayons et surface
    ri = Di / 2
    ro = Do / 2
    A_ext = math.pi * Do * L  # Surface externe totale

    # 1. Calcul initial de h_i (indépendant de T_wall)
    h_i = h_internal(m_dot, Di, T_coil_avg, L)

    # 2. Résistance thermique de la paroi
    R_wall = ro * math.log(ro / ri) / k_wall  # [m²·K/W]

    # 3. Estimation initiale de la température de paroi
    T_wall = 0.5 * (T_coil_avg + T_tank)

    # 4. Itération pour convergence
    for i in range(max_iter):
        # Calcul de h_o avec estimation courante de T_wall
        h_o = h_external(Do, T_wall, T_tank)

        # Calcul de U basé sur la surface externe
        U = 1 / (1/h_o + R_wall + (ro/ri)/h_i)

        # Calcul de la nouvelle T_wall via bilan d'énergie
        Q_dot = U * A_ext * (T_tank - T_coil_avg)  # Flux thermique total [W]
        T_wall_new = T_tank - Q_dot / (h_o * A_ext)

        # Vérifier convergence
        if abs(T_wall_new - T_wall) < tol:
            return U, Q_dot, T_wall_new, h_i, h_o

        # Relaxation pour stabilité
        T_wall = 0.7 * T_wall + 0.3 * T_wall_new

    # Si pas de convergence
    print(f"Attention : pas de convergence après {max_iter} itérations")
    return U, Q_dot, T_wall, h_i, h_o

# =============================================
# EXEMPLE D'UTILISATION
# =============================================
if __name__ == "__main__":
    # Paramètres du système
    params = {
        'Di': 0.009,          # Diamètre intérieur tube [m]
        'Do': 0.010,          # Diamètre extérieur tube [m]
        'k_wall': 15.0,       # Conductivité inox [W/m·K]
        'm_dot': 200/3600,    # Débit massique [kg/s]
        'T_coil_avg': 78.0,   # Température moyenne eau dans serpentin [°C]
        'T_tank': 92.0,       # Température eau dans la cuve [°C]
        'L': 10.0             # Longueur totale du serpentin [m]
    }

    # Calcul
    U, Q_dot, T_wall, h_i, h_o = calculate_U(params)

    # Résultats
    print("=============================================")
    print("RÉSULTATS DU CALCUL")
    print("=============================================")
    print(f"• Coefficient global U: {U:.1f} W/m²·K")
    print(f"• Flux thermique total: {Q_dot:.1f} W")
    print(f"• Température de paroi estimée: {T_wall:.1f} °C")
    print(f"• Coefficient interne h_i: {h_i:.1f} W/m²·K")
    print(f"• Coefficient externe h_o: {h_o:.1f} W/m²·K")
    print("=============================================")