import numpy as np
import math
from CoolProp.CoolProp import PropsSI

# Calculer échanges de chaleur et charges

# Valider avec https://www.pressure-drop.com/Online-Calculator/dp.php
# Exemple de résultat:
# Pressure Drop Online-Calculator
# Calculation output
# Flow medium: 	Cream 65°C / liquid
# Volume flow: 	0.25 m3/h
# Weight density: 	1020 kg/m3
# Dynamic Viscosity: 	9 mPa s
#
# Element of pipe: 	circular
# Dimensions of element: 	Diameter of pipe D: 9.5 mm
# Length of pipe L: 30 m
#
# Velocity of flow: 	0.98 m/s
# Reynolds number: 	1055
# Velocity of flow 2: 	-
# Reynolds number 2: 	-
# Flow: 	laminar
# Absolute roughness: 	0.1 mm
# Pipe friction number: 	0.06
# Resistance coefficient: 	191.6
# Resist.coeff.branching pipe: 	-
# Press.drop branch.pipe: 	-
# Pressure drop: 	937.92 mbar = 0.94 bar

# Pascals in a Bar, milli Bar...
PaInBar = 1e5
PaInmBar = 100
mmInM = 1000

pasto = "nefis"
debitsLH = [50,100,150,200,250,300]

FLUID_WATER = "Eau"
FLUID_MILK = "Lait"
FLUID_CREAM = "Crème"

# ======================
# PARAMÈTRES INITIAUX
# ======================
raccords = 18 # raccords rétrécisants de 9,5mm à 8,7mm
coudes = 8 # coudes à 90° (connexions aux échangeurs et à certains capteurs de température)

fluid = FLUID_MILK
diam_chauffe_ext = None # if not specified, will be calculated...
paroi_chauffe_defaut = 0.0005 # 0.5mm walls by default for coils giving the additional heat
temp_chauffe = 85
duree_maintien = 15 # secondes
longueur_maintien = None # if not specified, will be calculated...
convection_bain_marie = 500	# W / m ² °C. Estimation de transmission de chaleur par les courants convectifs naturels dans un bain marie
diam_calandre = None

a_plaque = None
parallele = 1 # Nb. d'échanges en parallèle (tuyaux multiples ou canaux multiples de l'échangeur à plaques)
diam_calandre_chauffe = None
debit_circulateur = None
diam_serpentin = 0.5 # Grand serpentin! Mettre la bonne valeur

diam_jonction = 0.0095
long_jonction = 2 #juste après la pompe et jusque avant la sortie, 8 sinon
obstacles = None

temp_nettoyage = 60
cleaning_fluid = FLUID_WATER #  Pour une solution de soude caustique à 1%, la différence de viscosité par rapport à l'eau est minime et souvent négligeable dans la pratique
v_nettoyage = 1.5 # m/s
debit_nettoyage = None

if pasto == "2022":
    debit = 0.150 / 3600  # m³/s (300 L/h)
    temp_froid_in = 4         # °C (lait entrant)
    temp_chaud_in = 85        # °C (lait sortant du pasteurisateur)

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.0095    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.0105    # m (diamètre extérieur tube central)
    diam_calandre = 0.014    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.28  # Diamètre de courbure du serpentin (m)

    longueur = 18 # hypothese de départ

    appoint = 2380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = 0.008
    longueur_chauffe = 28
    diam_maintien = 0.0095

elif pasto == "2023":
    debit = 0.2 / 3600  # m³/s (300 L/h)
    temp_froid_in = 4         # °C (lait entrant)
    temp_chaud_in = 85        # °C (lait sortant du pasteurisateur)
    temp_chauffe = 89

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.009    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.010    # m (diamètre extérieur tube central)
    diam_calandre = 0.013    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.35  # Diamètre de courbure du serpentin (m)

    longueur = 21 # hypothese de départ

    appoint = 2380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = 0.008
    diam_chauffe_ext = 0.009
    longueur_chauffe = 28
    diam_maintien = 0.0095
    longueur_maintien = 12

elif pasto == "upgrade2022":
    debitsLH = [70,
                144,
                178,
                210,
                250,
                282,
                312 ]
    debit = 0.250/ 3600  # m³/s (300 L/h)
    temp_froid_in = 20        # °C (lait entrant)
    temp_chaud_in = 85        # °C (lait sortant du pasteurisateur)
    temp_froid_in = 21        # °C (lait entrant)
    temp_chaud_in = 73.3        # °C (lait sortant du pasteurisateur)
    temp_chauffe = 89
    temp_chauffe = 75.3

    fluid = FLUID_WATER
    temp_froid_in = 21        # °C (lait entrant)
    temp_chaud_in = 23        # °C (lait sortant du pasteurisateur)
    temp_chauffe = 25

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.008#0.0095    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.010#0.0105    # m (diamètre extérieur tube central)
    diam_calandre = 0.014    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.28  # Diamètre de courbure du serpentin (m)

    longueur = 18 # 2 x 9 mètres

    appoint = 4380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = 0.009
    longueur_chauffe = 15
    diam_chauffe = (((0.0113**2 * 6) + (0.0095**2 *9))/15)**0.5   #0.0117
    diam_chauffe_ext = (((0.0127**2 * 6) + (0.0105**2 * 9))/15)**0.5 #0.0125
    longueur_chauffe = 15 #20

    diam_maintien = 0.0095
    longueur_maintien = 9

elif pasto == "12bonniers":
    debit = 0.250 / 3600  # m³/s (300 L/h)
    temp_froid_in = 4         # °C (lait entrant)
    temp_chaud_in = 85        # °C (lait sortant du pasteurisateur)

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.0095    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.0105    # m (diamètre extérieur tube central)
    diam_calandre = 0.014    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.28  # Diamètre de courbure du serpentin (m)

    longueur = 20 # hypothese de départ

    appoint = 4300 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = ((0.0113*6) + (0.0095*9))/15
    diam_chauffe_ext = ((0.0127*6) + (0.0105*9))/15
    longueur_chauffe = 15 #9+ (0.2 * np.pi * 17.5)
    temp_chauffe = 89

    diam_maintien = 0.0095

elif pasto == "ideal":
    #Deux Coolosus 2 pour l'échangeur (récupérateur)
    #Un Coolosus 1 du même fabricant (6m) pour la circulation d'eau chaude
    debit = 0.5 / 3600  # m³/s (300 L/h)
    temp_froid_in = 4         # °C (lait entrant)
    temp_chaud_in = 85        # °C (lait sortant du pasteurisateur)

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.0113    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.0127    # m (diamètre extérieur tube central)
    diam_calandre = 0.020    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.28  # Diamètre de courbure du serpentin (m)

    longueur = 18 # hypothese de départ

    appoint = 8000 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = 0.0113
    diam_chauffe_ext = 0.0127
    diam_calandre_chauffe = 0.020
    longueur_chauffe = 6
    debit_circulateur = debit * 2
    temp_chauffe = 89

    diam_jonction = 0.013
    diam_maintien = 0.013

elif pasto == "nefis-cata":
    fluid = FLUID_CREAM
    debit = 0.090 / 3600  # m³/s (300 L/h)
    temp_froid_in = 55         # °C (lait entrant)
    temp_chaud_in = 89        # °C (lait sortant du pasteurisateur)

    # Dimensions initiales (à optimiser)
    diam_int_tube = 0.009 #119    # m (diamètre intérieur tube central)
    diam_ext_tube = 0.0105 #27    # m (diamètre extérieur tube central)
    diam_calandre = 0.014#21    # m (diamètre intérieur du tube extérieur)

    diam_serpentin = 0.28  # Diamètre de courbure du serpentin (m)

    longueur = 18 # hypothese de départ

    appoint = 2400 #4380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_chauffe = 0.008
    longueur_chauffe = 28
    temp_chauffe = 92

    diam_maintien = 0.0095
    longueur_maintien = 10

elif pasto == "nefis":
    debitsLH = [37,
                74,
                125,
                161,
                189,
                216,
                268,
                313,
                342,
                370
                ]

    #Échangeur à plaques
    #Tube dans tube "Klarstein" avec circulation forcée
    fluid = FLUID_CREAM
    a_plaque = True
    debit = 0.250 / 3600  # m³/s (300 L/h)
    temp_froid_in = 55         # °C (lait entrant)
    temp_chaud_in = 89        # °C (lait sortant du pasteurisateur)

    # fluid = FLUID_WATER
    # temp_froid_in = 21        # °C (lait entrant)
    # temp_chaud_in = 23        # °C (lait sortant du pasteurisateur)
    # temp_chauffe = 25

    # Paramètres géométriques
    plaques_serie = 1
    parallele = 14    # nombre d'échanges en parallele (et non en série)
    longueur_plaque = 0.35  # m, longueur d'une plaque
    largeur_plaque = 0.11
    surface_plaque = longueur_plaque * largeur_plaque #0.043  # m², surface d'une plaque
    angle_chevron = 60  # degrés (typique)
    epaisseur_plaque = 0.0007 #m = 1.5mm d'épaisseur de plaque

    espacement = (0.100 / 29) - epaisseur_plaque #m = épaisseur totale / nombre de plaques #était 0.0021 #volume_plaque / surface_plaque
    # Calcul dérivés
    volume_plaque = surface_plaque * espacement #0.0002
    volume_canal = volume_plaque * plaques_serie * parallele # m³ (200 mL), volume d'une sectionl
    longueur = longueur_plaque * plaques_serie
    section_canal = espacement * largeur_plaque #volume_plaque / longueur_plaque

    appoint = 4380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_serpentin = 0.28 # pour l'appoint
    diam_chauffe = 0.0113 #0.0095 #0.0113
    diam_ext_chauffe = 0.0127 #0.0127    # m (diamètre extérieur tube central)
    #diam_calandre_chauffe = 0.014 #0.020    # m (diamètre intérieur du tube extérieur)
    #debit_circulateur = debit * 3
    longueur_chauffe = 15.2
    temp_chauffe = 93

    diam_maintien = 0.013
    longueur_maintien = 9

    diam_calandre_chauffe = None
    debit_circulateur = None

    diam_jonction = 0.013

elif pasto == "nefis-eau":

    #Échangeur à plaques
    a_plaque = True
    debit = 0.268 / 3600  # m³/s (300 L/h)
    fluid = FLUID_WATER
    temp_froid_in = 21
    temp_chaud_in = 23

    # Paramètres géométriques
    plaques_serie = 1
    parallele = 14    # nombre d'échanges en parallele (et non en série)
    longueur_plaque = 0.35  # m, longueur d'une plaque
    largeur_plaque = 0.11
    surface_plaque = longueur_plaque * largeur_plaque #0.043  # m², surface d'une plaque
    angle_chevron = 60  # degrés (typique)
    epaisseur_plaque = 0.0007 #m = 1.5mm d'épaisseur de plaque

    espacement = (0.100 / 29) - epaisseur_plaque #m = épaisseur totale / nombre de plaques #était 0.0021 #volume_plaque / surface_plaque
    # Calcul dérivés
    volume_plaque = surface_plaque * espacement #0.0002
    volume_canal = volume_plaque * plaques_serie * parallele # m³ (200 mL), volume d'une sectionl
    longueur = longueur_plaque * plaques_serie
    section_canal = espacement * largeur_plaque #volume_plaque / longueur_plaque

    appoint = 4380 # W en supposant 500W de pertes (2880W = 14A + 1920W = 8A)4800W d'électricité à 20A)

    diam_serpentin = 0.28 # pour l'appoint
    longueur_chauffe = 9.4
    temp_chauffe = 93
    #TEST
    diam_chauffe = ((0.0113*6) + (0.0095*9))/15   #0.0117
    diam_chauffe_ext = ((0.0127*6) + (0.0105*9))/15 #0.0125
    longueur_chauffe = 15
    diam_calandre_chauffe = None
    debit_circulateur = None

    diam_jonction = 0.0095
    diam_maintien = 0.0095
    longueur_maintien = 9

if a_plaque is None:
    epaisseur_plaque = (diam_ext_tube - diam_int_tube) / 2 # wall
if diam_chauffe_ext is None:
    diam_chauffe_ext = diam_chauffe + paroi_chauffe_defaut # 0.5mm walls by default for coils giving the additional heat
if debit_circulateur is None:
    debit_circulateur = debit
if longueur_maintien is None:
    longueur_maintien = duree_maintien * debit / ((diam_maintien * diam_maintien / 4) * np.pi)

if obstacles is None: # Définir les obstacles
        obstacles = [
            ('coudes', 0.4, diam_jonction, coudes),      # 6 coudes
            ('retrecissement', diam_jonction, diam_jonction-0.0008, raccords),  # 18 rétrécissements 9.5→8.7mm
            ('elargissement', diam_jonction-0.0008, diam_jonction, raccords)    # 18 élargissements 8.7→9.5mm
        ]

# ======================
# FONCTIONS DE CALCUL
# ======================
def proprietes_fluide(fluid: str, temp_C: float) -> dict:
    """Retourne les propriétés physiques du lait à température donnée"""
    temp_K = temp_C + 273.15
    if fluid == FLUID_WATER:
        return {
            # Get the density of Water at T = temp (in Kelvin) and P = 101325 Pa (ground air pressure)
            'rho': PropsSI('D', 'T', temp_K, 'P', 101325, 'Water'), # kg/m³
            # thermal capacity in Joule
            'cp': PropsSI('C', 'T', temp_K, 'P', 101325, 'Water'), # J / kg / K
            # thermal conductivity
            'k': PropsSI('L', 'T', temp_K, 'P', 101325, 'Water'), # W/m/K
            # Viscosity
            # 20°C, water's dynamic viscosity is about 1.002 mPa·s, which is equivalent to 0.001 Pa.s or 1 centipoise
            'mu': PropsSI('V', 'T', temp_K, 'P', 101325, 'Water'), # Pa . s
            # PRANDTL number
            'pr': PropsSI('PRANDTL', 'T', temp_K, 'P', 101325, 'Water')
        }
        # 'rho': 1000 - 0.013 * temp_c**2,  # kg/m³
        # 'cp': 4180 + 1.5 * temp_c,        # J/kg·K
        # 'k': 0.6 + 0.0017 * temp_c,       # W/m·K
        # 'mu': 0.001 * np.exp(-0.03*temp_c),# Pa·s (viscosité dynamique)
        # 'pr': 6.5 * np.exp(-0.05*temp_c)  # Nombre de Prandtl
        # Plus précisément:
        #'rho': 997.0 * (1 - (temp_c - 4)**2 * 2e-5),  # kg/m³
        #'cp': 4181.7 - 2.72 * temp_c + 0.013 * temp_c**2,  # J/kg·K
        #'k': 0.606 * (1 + 0.001 * temp_c),  # W/m·K
        #'mu': 0.001 * np.exp(-1.94 - 4.80/(temp_c + 42.5))  # Pa·s
        # ou encore en Kelvin:
        # T = temp_c + 273.15
        # # Masse volumique (kg/m³)
        # rho = 999.9 - 0.0226*T + 2.87e-5*T**2 - 1.13e-8*T**3
        # # Chaleur spécifique (J/kg·K)
        # cp = 4217.4 - 3.720*T + 0.0152*T**2 - 2.0e-5*T**3
        # # Conductivité thermique (W/m·K)
        # k = 0.558 + 0.00216*T - 1.03e-5*T**2
        # # Viscosité dynamique (Pa·s) - Formule de Vogel
        # mu = 0.001 * np.exp(-3.7188 + 578.919/(T-137.546))
        # pr = (cp * mu) / k
    #elif fluid == "milk": # 4%, 8,95% solids not fat
    #    • Masse volumique du lait (ρ) ≈ 1030 kg/m³ (légèrement plus dense que l'eau)
    # Viscosité dynamique du lait (μ) ≈ 2,0 × 10⁻³ Pa·s (pour du lait à 4% de matière grasse à 20°C)
    elif fluid == FLUID_MILK:  # cru: 3.8% MG
        MG = 3.8
        # Masse volumique [kg/m³]
        rho = 1032.5 - 0.275 * temp_C - 0.0015 * temp_C**2
        #rho = 1034.5 + 0.3*temp_C - 0.03*temp_C**2 - 0.7*MG - (0.01 * MG**2)

        # Viscosité dynamique [mPa·s]
        #mu = 2.06e-3 * np.exp(3425 / temp_K)
        #mu = np.exp( -8.9 + (0.1 *3.8) + (2721.5 / temp_K))
        mu = 226.9503858568 * np.exp(-0.0162557770797138 * temp_K) / 1000 #Pa.s

        # Chaleur spécifique [J/(kg·K)]
        cp = 3650 + 6.5 * temp_C - 0.02 * temp_C**2

        # Conductivité thermique [W/(m·K)]
        k = 0.535 + 0.0014 * temp_C

        return {
            'rho': rho,
            'mu': mu,  # Converti en Pa·s
            'cp': cp,
            'k': k,
            'pr': (cp * mu) / k
        }
    elif fluid == FLUID_CREAM: # 37%
        MG = 37 #%
        #rho0 = 1048.2 - 0.488*temp_C + 0.0012*temp_C**2
        rho = 1034.5 + 0.3*temp_C - 0.03*temp_C**2 - 0.7*MG - (0.01 * MG**2)
        # Viscosité dynamique (Pa.s)
        #mu0 = 4.85e-7 * np.exp(3279.4 / temp_K)# 20°C=15-115 mPas
        mu = np.exp( -8.9 + (0.1 * MG) + (2721.5 / temp_K)) / 1000

        # Chaleur spécifique (J/kg.K)
        cp = 3180 + 8.5*temp_C
        # Conductivité thermique (W/m.K)
        k = 0.382 + 0.0012*temp_C - 5e-6*temp_C**2
        return{
            # # Masse volumique (kg/m³)
            'rho': rho,
            # # Chaleur spécifique (J/kg·K)
            'cp': cp,
            # # Conductivité thermique (W/m·K)
            'k': k,
            # # Viscosité dynamique (Pa·s) - Formule de Vogel
            'mu': mu,
            'pr': (cp * mu) / k
        }

print (f"Propriétés de {fluid}")
print (" \tMasse\tChaleur\tConduc\tViscoD\tPrank")
print ("°C\tkg/m³\tJ/kg.K\tW/m.K\tmPa.s\t")
for t in [20,40,60,80]:
    props = proprietes_fluide(fluid, t)
    print (f"{t}\t{props['rho']:.1f}\t{props['cp']:.1f}\t{props['k']:.3f}\t{props['mu']*1000:.3f}\t{props['pr']:.1f}")

# =============================================
# DONNÉES PHYSIQUES (CONSTANTES)
# =============================================
k_inox = 15.0  # Conductivité thermique de l'inox (W/m·K)

# ======================
# BILAN THERMIQUE
# ======================

def calcul_lmtd(Tci, Tco, Thi, Tho):
    """Calcule le LMTD pour échangeur à contre-courant"""
    # Source: http://www.heatgroup.eesc.usp.br/tools/5210.php#results
    delta_T1 = Thi - Tco
    delta_T2 = Tho - Tci
    return (delta_T1 - delta_T2) / np.log(delta_T1 / delta_T2) if delta_T1 != delta_T2 else delta_T1

def nombre_dean(Re: float, d_h: float, D_curve: float) -> float:
    """Calcule le nombre de Dean pour les effets de courbure"""
    # Validé par les interprétations de Dean, W.R. (1927). "Note on the motion of fluid in a curved pipe". Philosophical Magazine, 4(20), 208-223.*
    # mais sa lecture est difficile:
    # https://www.cambridge.org/core/journals/mathematika/article/note-on-the-motion-of-fluid-in-a-curved-pipe/C53740B97524473762A91F32F38E8781
    dean = Re * ((d_h / D_curve)**0.5)
    print (f"Re={Re:.1f}, rac d/D={(d_h/D_curve)**0.5:.5f}, Dean={dean:.1f}")
    return dean

def correction_serpentin_nu(Nu_straight: float, De: float, section_type: str) -> float:
    """Applique la correction de Nusselt pour les serpentins"""
    # Référence de DeepSeek:
    # Schmidt, E.F. (1967). "Wärmeübergang und Druckverlust in Rohrschlangen". Chemie Ingenieur Technik, 39(13), 781-789
    # if section_type == 'tube':
    #     # Corrélation de Schmidt pour tubes courbés
    #     return Nu_straight * (1 + 0.14 * De**0.7)
    # else:
    #     # Corrélation approximative pour annulaire
    #     return Nu_straight * (1 + 0.075 * De**0.6)
    """Correction réaliste avec saturation pour De élevés"""
    # Limite De à 200 pour éviter des corrections irréalistes
    De_clamped = De#min(De, 200)

    if section_type == 'tube' or section_type == 'maintain':
        # Corrélation de Schmidt modifiée
        return Nu_straight * (1 + 0.8 * (1 - np.exp(-0.05 * De_clamped)))
    else:
        # Corrélation pour annulaire avec saturation
        return Nu_straight * (1 + 0.6 * (1 - np.exp(-0.03 * De_clamped)))

# def correction_serpentin_f(f_straight: float, De: float) -> float:
#     """Applique la correction du facteur de friction"""
#     # return f_straight * (0.37 * De**0.36)
#     """Correction réaliste pour le facteur de friction"""
#     De_clamped = De# min(De, 200)
#     # Référence de DeepSeek:
#     # Ito, H. (1959). "Friction Factors for Turbulent Flow in Curved Pipes". Journal of Basic Engineering, 81(2), 123-134.
#     return f_straight * (1 + 0.3 * (1 - np.exp(-0.02 * De_clamped)))

# =============================================
# FONCTIONS DE CALCUL
# =============================================
def resistance_paroi(diam_int_tube: float, diam_ext_tube: float) -> float:
    """Calcule la résistance thermique de la paroi (m²K/W)"""
    # Calcul de l'épaisseur et de la surface moyenne
    epaisseur = (diam_ext_tube * np.log(diam_ext_tube/diam_int_tube)) / 2

    # Résistance thermique = épaisseur / conductivité
    return epaisseur / k_inox #  * np.log(diam_ext_tube/diam_int_tube)

# def calcul_pertes_charge_serpentin(L: float, d_h: float, f: float, rho: float, v: float) -> float:
#     """Calcule les pertes de charge avec correction serpentin"""
#     # Pertes de charge régulières
#     dp_linear = f * (L / d_h) * (rho * v**2 / 2)
#
#     # Pertes supplémentaires dues aux coudes (1 coude par mètre)
#     n_coude = L  # Approximation
#     k_coude = 0.5 * (d_h / diam_serpentin)**0.5  # Coefficient de perte localisée
#     dp_coude = n_coude * k_coude * (rho * v**2 / 2)
#
#     return dp_linear + dp_coude

def coefficient_convection_serpentin(fluid: str, parallele: int,
                                     diam_int_tube: float ,diam_ext_tube: float, diam_calandre: float, longueur: float,
                                     diam_serpentin: float, debit: float, temp: float, section_type: str) -> tuple:

    props = proprietes_fluide(fluid, temp)

    """Coefficient de convection avec correction serpentin"""
    # Calcul section Area
    if section_type == 'tube' or section_type == "maintain":
        A = np.pi * (diam_int_tube/2)**2
        d_h = diam_int_tube
    else:
        A = np.pi * (diam_calandre**2 - diam_ext_tube**2) / 4.0
        d_h = diam_calandre - diam_ext_tube

    v = (debit/parallele) / A
    print(f"Débit:{debit:.6f} m³/s / {parallele}, Vitesse: {v:.1f} m/s")

    Re = props['rho'] * v * d_h / props['mu']
    #print(f"Temp={temp:.1f}, Rho={props['rho']:.1f}, Diam.Hydro={d_h:.3f}, µ={props['mu']:.6f}, Reynold: {Re:.1f}")

    # Calcul de base (tube droit) basé sur:
    # Gnielinski, V. New equations for heat and mass transfer in turbulent pipe and channel flow. Int Chem Eng 16:359-367. (1976)
    # et testable sur http://www.heatgroup.eesc.usp.br/tools/1114.php
    # Référence très complète: https://slideplayer.com/slide/12099829/

    # Re = (ρ * v * D) / μ
    #
    # if Re < 2000:
    #     f = 64 / Re  # Régime laminaire
    # else:
    #     f = 0,035   # Régime turbulent
    #
    # ΔP = f * (L / D) * (ρ * v² / 2) + ΣK * (ρ * v² / 2)  # Pertes singulières

    Nu2=0
    diam_ratio = d_h/diam_serpentin
    De = nombre_dean(Re, d_h, diam_serpentin)
    Nu_laminaire0 = 3.66 * (1+ (0.049*(De**0.75)))
    Nu_turbulent0 = 0.023*(Re**0.8)*(props['pr']**0.4)*(1+(0.1*(min(De,500)**0.5)))

    Nu_laminaire = 0.7 * (Re**0.43)*(props['pr']**(1/6)) * (diam_ratio**0.07)
    Nu_turbulent = 0.00619 * (Re**0.92) * (props['pr']**0.4) * (1 + 3.46*diam_ratio)
    if Re > 2000:
        f = (0.79 * np.log(Re) - 1.64)**-2
        Nu = (f/8) * (Re - 1000) * props['pr'] / (1 + 12.7 * ((f/8)**0.5) * (props['pr']**(2/3) - 1))
        #Nu2 = 0.023 * Re**0.8 * props['pr']**0.4
    else:
        #Nu2 = 3.66 + (0.0668 * (diam_int_tube/longueur) * Re * props['pr']) \
        #      / (1 + 0.04 * ((diam_int_tube/longueur) * Re * props['pr'])**(2/3))
        if section_type == 'tube' or section_type == "maintain":
            Nu = 3.66 + 0.19 * (Re * props['pr'] * d_h/longueur)**0.8
        else:
            Nu = 4.86 + 0.27 * (Re * props['pr'] * d_h/longueur)**0.5
        f = 64 / Re if Re > 0 else 0
    # Correction pour serpentin
    #print (f"\n{section_type}: Reynold={Re:.1f}, Dean for dh={d_h*1000:.1f}mm : {De:.1f}")
    Nu_curved = correction_serpentin_nu(Nu, De, section_type)
    Nu_Radwan = 0.036*(min(De,500)**0.85)*(props['pr']**0.4)*(1+(0.2*min(De,500)**0.5))
    print (f"Nu={Nu:.3f}, Nu2={Nu2:.3f}, Nu curved={Nu_curved:.3f}, Radwan Nu={Nu_Radwan:.3f}, Laminaire={Nu_laminaire:.1f}, Turbulent={Nu_turbulent:.1f}")
    print (f"NuT={Nu_turbulent0:.3f}, NuL={Nu_laminaire0:.3f}, Radwan Nu={Nu_Radwan:.3f}, Laminaire={Nu_laminaire:.1f}, Turbulent={Nu_turbulent:.1f}")    # if section_type == 'tube':
    #      print(f"Direct Nu={0.021*(De**0.8)*(props['pr']**0.4):.3f}")
    # else:
    #      print(f"Direct Nu={0.0245*(De**0.77)*(props['pr']**0.4):.3f}")
    # print(f"Empirical Nu={0.023*(Re**0.8)*(1.0+(0.1*(De**0.5)))*(props['pr']**0.4):.3f}")

    # f_curved = correction_serpentin_f(f, De)
    if section_type == 'tube' or section_type == "maintain":
        # gamma = d_h / diam_serpentin
        # f_turbulent = 0.3164 * (Re ** -0.25) + 0.012 * gamma**0.5

        # # Apply roughness correction if applicable (Das's correlation Eq. 31)
        if section_type == 'tube':
             roughness = 4.5e-5  # m (typical for stainless steel)
        else: #elif material == 'LDPE':
             roughness = 1e-6    # m (smooth plastic)
        # if roughness > 0:
        #     F_roughness = 4*17.5782 * (Re ** -0.3137) * (gamma ** 0.3621) * ((roughness / d_h) ** 0.6885)
        #     f_turbulent += F_roughness

        kA = -2 * math.log10(roughness/d_h / 3.7 + 12/Re)
        kB = -2 * math.log10(roughness/d_h / 3.7 + 2.51*kA/Re)
        kC = -2 * math.log10(roughness/d_h / 3.7 + 2.51*kB/Re)
        f_turbulent = (kA - (kB - kA)**2 / (kC - 2*kB + kA))**(-2) * (1+(0.033*(min(De,500)**0.5)))
        #f_turbulent = (0.3164 / (Re**0.25)) * min(1 + (0.01 * De), 3.0 )  #0.01 au lieu de 0.033    +(0.03*(De**-0.3))
    else:
        #f_turbulent = 0.079*(Re**-0.25)*(1+(0.1*(De**0.5)))
        k = diam_ext_tube / diam_calandre
        denom_corr = (1 - (k**2) - ((k**4) / 2))
        f_turbulent = (0.0791 * (Re**-0.25)) * ((1 - k**2) / denom_corr) * (1+(0.033*(min(De,500)**0.5)))
        # stainless_coef = 0.015
        # A = (2.457 * np.log(1/((7/Re)**0.9 + 0.27*stainless_coef)))**16
        # B = (37530/Re)**16
        # f_turbulent = 8 * ((8/Re)**12 + 1/(A+B)**1.5)**(1/12)
    f_laminaire = (64/Re) * (1+(0.1*(min(De,500)**0.5)))

    if Re > 5000: # and De > 500:
        Nu_final = Nu_Radwan
        f_interpol = f_turbulent
    elif Re >= 2000:
        Nu_final = (Nu_laminaire**3+Nu_Radwan**3)**(1/3)
        X = (Re-2000)/3000
        f_interpol = (f_laminaire*(1.0-X)) + (f_turbulent*X)
    else:
        Nu_final = Nu_laminaire
        f_interpol = f_laminaire
    #print (f"Nu retenu: {Nu_final:.1f}")
    #print(f"Calculs de \"f\", Actuel={f:.6f}, Curved={f_curved:.6f}, Turbulent={f_turbulent:.6f}, Laminaire={f_laminaire:.6f}, Retenu={f_interpol:.6f}")

    # Calcul des pertes de charge côté chaud
    #delta_P = calcul_pertes_charge_serpentin(longueur, d_h, f_interpol, props['rho'], v)
    delta_P2 = None
    delta_P = f_interpol * (longueur / d_h) * (props['rho'] * v**2 / 2)

    h = Nu_final * props['k'] / d_h
    return h, f_interpol, Re, De, v, delta_P, delta_P2

def coefficient_convection_plaque(fluid: str, section_canal: float, longueur: float, angle_chevron: int, debit: float,
                                  temp: float, parallele: float, diam_jonction: float) -> tuple:

    props = proprietes_fluide(fluid, temp)
    # https://achp.sourceforge.net/ACHPComponents/PlateHeatExchanger.html
    #d_h = 2 * espacement
    d_h = 4 * section_canal / ((largeur_plaque+espacement)*2)
    print(f'Hydraulic diameter={d_h*mmInM:.1f}mm')
    v = (debit/parallele) / section_canal
    print(f'v={v*mmInM:.1f} mm/sec')

    Re = props['rho'] * v * d_h / props['mu']
    print(f'Reynold={Re:.1f}')

    # Calcul de base (tube droit) basé sur:
    # Gnielinski, V. New equations for heat and mass transfer in turbulent pipe and channel flow. Int Chem Eng 16:359-367. (1976)
    # et testable sur http://www.heatgroup.eesc.usp.br/tools/1114.php
    # Référence très complète: https://slideplayer.com/slide/12099829/
    if Re > 3000:
        Nu = 0.28 * (Re**0.65) * (props['pr']**0.36) * (1 + 0.0015*angle_chevron**1.2)
    else: # creme et autres liquides visqueux
        Nu = 0.22 * (Re**0.6) * (props['pr']**0.4) * (1 + 0.0015*angle_chevron**1.2)

    h = Nu * props['k'] / d_h

    A = np.pi * (diam_jonction/2)**2  # Section étroite
    vj = debit / A

    C = 0.26
    g = 0.25
    f = C * (Re ** g)  # corrélation empirique (Eq. 3.38)
    #delta_P_canal = (4 * f * longueur / d_h) * 0.5 * props['rho'] * v**2  # perte de charge (Eq. type Darcy/Fanning)

    # # Facteur friction (corrélation Martin)
    # f = (1.5 / Re**0.5) * (1 + 0.002*angle_chevron**1.5)
    # print(f"f={f:.1f}")
    # #else: # Low Re, short channel: https://research.library.mun.ca/15388/1/thesis.pdf
    # if Re < 2300:
    #     f = (12 / Re)\
    #         * (((2/(math.cos((math.radians(angle_chevron)))**1.73))**2)
    #            +((1.33/((longueur/d_h)**0.5))*(Re**(0.0551*(angle_chevron**0.675))))**2)**0.5
    #     print(f"f low Re={f:.1f}")
    # Pertes de charge totales
    delta_P_canal = 100 * f * (longueur/d_h) * (props['rho'] * v**2)
    print(f"canal={delta_P_canal:.1f}Pa")
    #delta_P_entree_sortie = 1.4 * (props['rho'] * v**2 / 2)

    beta = diam_jonction / (diam_jonction*parallele)
    K = 1.5 * (1 - beta**2)
    delta_P_entree_sortie = K * (props['rho'] * vj**2) / 2
    print(f"inlet/outlet={delta_P_entree_sortie:.1f}Pa")
    delta_P = delta_P_canal + delta_P_entree_sortie

    return h, f, Re, v, delta_P

def coefficient_convection_trempe (fluid, Di: float, Do: float, longueur: float, diam_serpentin: float,
                                   temp_tube: float, temp_chauffe: float):

    props = proprietes_fluide(fluid, temp_chauffe)
    beta = 0.0002 * np.exp(0.005 * temp_chauffe)
    nu = props['mu'] / props['rho']  # Viscosité cinématique [m²/s]
    alpha = props['k'] / (props['rho'] * props['cp'])  # Diffusivité thermique [m²/s]
    T_wall = temp_tube
    T_mean = 0.5 * (T_wall + temp_chauffe)
    delta_T = abs(T_wall - temp_chauffe)

    # Éviter delta_T = 0 pour les calculs
    if delta_T < 0.1:
        delta_T = 0.1

    Gr = (9.81 * beta * delta_T * Do**3) / (nu**2)
    Ra = Gr * props['pr']

    # Corrélation de Churchill-Chu
    term = (1 + (0.559 / props['pr'])**(9/16))**(8/27)
    Nu = (0.60 + 0.387 * Ra**(1/6) / term)**2
    return Nu * props['k'] / Do

def pertes_singulieres(Q, obstacles, rho):
    """
    Calcule les pertes singulières dynamiques
    Args:
        Q: Débit (m³/s)
        obstacles: Liste de tuples [('type', geo_params), ...]
                   ex: [('coudes', 0.4, D_ref), ('retrecissement', D1, D2)]
    Returns:
        ΔP_sing (bar)
    """
    delta_P = 0
    for obs in obstacles:
        type_obs = obs[0]

        if type_obs == 'coudes':
            K, D_ref, n = obs[1], obs[2], obs[3]
            A = np.pi * (D_ref/2)**2
            v = Q / A
            delta = n * K * (rho * v**2) / 2
            print (f"{type_obs}: {n:.0f}x {delta} Pa")
            delta_P += delta


        elif type_obs == 'retrecissement':
            D1, D2, n = obs[1], obs[2], obs[3]
            beta = min(D1, D2) / max(D1, D2)
            K = 0.5 * (1 - beta**2)
            A = np.pi * (min(D1, D2)/2)**2  # Section étroite
            v = Q / A
            delta = n * K * (rho * v**2) / 2
            print (f"{type_obs}: {n:.0f}x {delta} Pa")
            delta_P += delta

        elif type_obs == 'elargissement':
            D1, D2, n = obs[1], obs[2], obs[3]
            beta = min(D1, D2) / max(D1, D2)
            K = (1 - beta**2)**2
            A = np.pi * (min(D1, D2)/2)**2  # Section étroite
            v = Q / A
            delta = n * K * (rho * v**2) / 2
            print (f"{type_obs}: {n:.0f}x {delta} Pa")
            delta_P += delta

    print(f"Pertes singulières totales: {delta_P / PaInBar :.4f} bar")
    return delta_P  # Conversion Pa → bar

gravitational_constant = 6.6743e-11 # m3 kg-1 s-2

def charge_ligne_droite(fluid, debit, longueurs,temperatures, diametres):
    total_load = 0.0
    for i in range(len(longueurs)):
        long = longueurs[i]
        temp = temperatures[i]
        diam = diametres[i]
        props = proprietes_fluide(fluid, temp)
        A = np.pi * (diam/2)**2
        v = debit / A
        Rey = diam * v * props['rho'] / props['mu']
        f = 16 / Rey
        #f = 0.316 * (diam*v*props['rho']/props['mu'])**-0.25
        #fanning = 32 * f * debit**2 / (2 * gravitational_constant * diam**5 )
        #load = longueur * fanning
        load = 2 * f * long / diam * props['rho'] * v**2
        print(f"{temp:.1f}: {long:.1f}m = {doc_P(load,None)}")
        total_load += load
    return total_load

def doc_P(P1: float, P2: float) -> str:
    if P2 is None:
        return f"{P1 / PaInmBar :.1f}mbar"
    return f"{P1 / PaInmBar :.1f}mbar {((P2 - P1) / P1) * 100.0:.1f}%"

# =============================================
# EXEMPLE D'UTILISATION
# =============================================
# Paramètres géométriques


# =============================================
print("\n" + "="*50)
print(f"   Projet {pasto}")
print("="*50)
if a_plaque is None:
    print(f"Configuration: Tube dans tube (inox)")
    print(f"- Épaisseur paroi: {epaisseur_plaque*mmInM:.2f} mm")
    print(f"- Diamètre tube: {diam_int_tube*mmInM:.1f} mm int / {diam_ext_tube*mmInM:.1f} mm ext")
    print(f"  Longueur: {longueur*1000:.1f}mm, Volume: {longueur*diam_int_tube*diam_int_tube*math.pi*1000/4:.1f} L")
    print(f"- Diamètre calandre: {diam_calandre*mmInM:.1f} mm")
    print(f"  Longueur: {longueur*1000:.1f}mm, Volume: {longueur*(diam_calandre**2-diam_ext_tube**2)*math.pi*1000/4:.1f} L")
else:
    print(f"Configuration: Echangeur à 2x{plaques_serie*parallele}+1 plaques (inox)")
    print(f"- Plaques: {plaques_serie} séries de {parallele} en parallèle")
    print(f"- Épaisseur paroi: {epaisseur_plaque*mmInM:.2f} mm"
          f", {espacement*mmInM:.1f}mm d'espacement entre les plaques."
          f", Épaisseur totale: {((espacement+epaisseur_plaque)*(plaques_serie*parallele*2+1))*mmInM:.1f}mm")
    print(f"- Dimension: {longueur_plaque:.3f}m de longueur x {largeur_plaque:.3f}m de largeur,"
          f" surface={surface_plaque:.6f}m², volume={volume_plaque*1000:.3f}L,"
          f" côté complet={volume_canal*1000:.3f}L")
# !!! Résistances thermiques de la chauffe d'appoint
R_paroi_chauffe = resistance_paroi(diam_chauffe, diam_chauffe_ext)

results = []
for debit in debitsLH:
    debit = debit / (3600 * 1000)

    fluid_target = proprietes_fluide(fluid, temp_chaud_in)
    cp_fluid = fluid_target['cp']

    m_fluid = debit * fluid_target['rho']  # kg/s
    v_min = 1.5 # m/s Trouver la vitesse minimale dans le pasteurisateur
    temp_froid_out_target = temp_chaud_in - (appoint / (cp_fluid * m_fluid))
    if temp_froid_out_target < temp_froid_in:
        temp_froid_out_target = temp_froid_in

    # Températures moyennes
    tm_froid = (temp_froid_in + temp_froid_out_target) / 2

    # Flux thermique transféré (W)
    eau_mitigee = proprietes_fluide(fluid, tm_froid)
    m_fluid = debit * eau_mitigee['rho']  # kg/s
    cp_fluid = eau_mitigee['cp']
    q_transfert = m_fluid * cp_fluid * (temp_froid_out_target - temp_froid_in)

    temp_chaud_out = temp_chaud_in - (temp_froid_out_target - temp_froid_in )
    tm_chaud = (temp_chaud_in + temp_chaud_out) / 2

    lmtd = calcul_lmtd(temp_froid_in, temp_froid_out_target, temp_chaud_in, temp_chaud_out)
    # Températures moyennes
    temp_tube = (temp_chaud_in+temp_chaud_out)/2          # °C (eau chaude dans le tube)
    temp_calandre = (temp_froid_in+temp_froid_out_target)/2      # °C (eau froide dans la calandre)

    if a_plaque is None:
        # """Calcule le coefficient global de transfert thermique U (W/m²K)"""
        # Coefficients de convectionh, f_curved, Re, De, v
        h_tube, f_tube, rey_tube, dean_tube, v_tube, delta_P_tube, P2_tube = \
            coefficient_convection_serpentin(fluid, parallele,
                                             diam_int_tube, diam_ext_tube, diam_calandre, longueur, diam_serpentin, \
                                             debit, temp_tube, "tube")
        if v_tube < v_min:
            v_min = v_tube
        h_annulaire, f_annulaire, rey_annulaire, dean_annulaire, v_annulaire, delta_P_annulaire, P2_annulaire = \
            coefficient_convection_serpentin(fluid, parallele,
                                             diam_int_tube, diam_ext_tube, diam_calandre, longueur, diam_serpentin, \
                                             debit, temp_calandre, "shell")
        if v_annulaire < v_min:
            v_min = v_annulaire
        #print(f"Convection eau/tube int={h_tube:.3f}, ext={h_annulaire:.3f}")

        # !!! Résistances thermiques
        R_conv_tube = 1 / h_tube  # * np.pi * diam_int_tube_tube)
        R_conv_cal = 1 / h_annulaire # * np.pi * diam_ext_tube_tube)
        R_paroi = resistance_paroi(diam_int_tube, diam_ext_tube)
        #print(f"Résistances, int={R_conv_tube:.6f}, ext={R_conv_cal:.6f}, paroi={R_paroi:.6f}")

        A_ratio = diam_ext_tube / diam_int_tube  # Ratio des surfaces
        # Résistance totale et coefficient U
        R_total = (A_ratio*R_conv_tube) + R_paroi + R_conv_cal

        U_possible = 1 / R_total
        # ======================
        # SURFACE ET LONGUEUR
        # ======================
        surface_actuelle = (np.pi * diam_ext_tube * longueur)
        surface_echange = q_transfert / (U_possible * lmtd)
        longueur2 = surface_echange / (np.pi * diam_ext_tube)

    else:
        surface_actuelle = 2 * surface_plaque * parallele * plaques_serie
        # """Calcule le coefficient global de transfert thermique U (W/m²K)"""
        # Coefficients de convectionh, f_curved, Re, De, v
        h_tube, f_tube, rey_tube, v_tube, delta_P_tube = \
            coefficient_convection_plaque(fluid, section_canal, longueur, angle_chevron, debit, temp_tube, parallele, diam_jonction)
        delta_P_annulaire = delta_P_tube
        P2_tube = None
        P2_annulaire = None
        if v_tube < v_min:
            v_min = v_tube
        # !!! Résistances thermiques
        R_conv_tube = 1 / h_tube # * np.pi * diam_int_tube_tube)
        R_paroi = epaisseur_plaque / k_inox
        print(f"Résistances, chaque interface={R_conv_tube:.6f}, paroi={R_paroi:.6f}")

        # Résistance totale et coefficient U
        R_total = (2*R_conv_tube) + R_paroi

        U_possible = 1 / R_total
        # ======================
        # SURFACE ET LONGUEUR
        # ======================
        surface_echange = q_transfert / (U_possible * lmtd)
        longueur2 = surface_echange / largeur_plaque / parallele / 2

    h_chauffe, f_chauffe, rey_chauffe, dean_chauffe, v_chauffe, delta_P_chauffe, P2_chauffe = \
        coefficient_convection_serpentin(fluid, 1, diam_chauffe, diam_chauffe_ext,
                                         diam_calandre_chauffe, longueur_chauffe, diam_serpentin, \
                                         debit, (temp_chaud_in+temp_froid_out_target) / 2, "tube")
    if diam_calandre_chauffe is None:
        h_chauffe_annulaire = coefficient_convection_trempe(fluid, diam_chauffe, diam_chauffe_ext, longueur_chauffe, diam_serpentin, \
                                                             (temp_chaud_in+temp_froid_out_target)/2, temp_chauffe)
        delta_P_chauffe_annulaire = None
        P2_chauffe_annulaire = None
    else:
        h_chauffe_annulaire, f_chauffe_annulaire, rey_chauffe_annulaire, dean_chauffe_annulaire, v_chauffe_annulaire,\
        delta_P_chauffe_annulaire, P2_chauffe_annulaire = \
            coefficient_convection_serpentin(fluid, 1, diam_chauffe, diam_chauffe_ext,
                                             diam_calandre_chauffe, longueur_chauffe, diam_serpentin, \
                                             debit_circulateur, temp_chauffe, "shell")
    R_conv_chauffe = 1 / h_chauffe # * np.pi * diam_int_tube_tube)
    R_conv_chauffe_annulaire = 1 / h_chauffe_annulaire # * np.pi * diam_int_tube_tube)
    #print(f"Résistances appoint, int={R_conv_chauffe:.6f}, ext={1/convection_bain_marie:.6f}, paroi={R_paroi_chauffe:.6f}")
    # Résistance totale et coefficient U
    R_total_chauffe = R_conv_chauffe + R_paroi + R_conv_chauffe_annulaire
    U_possible_chauffe = 1 / R_total_chauffe

    h_maintien, f_maintien, rey_maintien, dean_maintien, v_maintien, delta_P_maintien, P2_maintien = \
        coefficient_convection_serpentin(fluid, 1, diam_maintien, None, None, longueur_maintien, diam_serpentin, \
                                         debit, temp_chaud_in, "maintain")
    if v_maintien < v_min:
        v_min = v_maintien

    h_jonction, f_jonction, rey_jonction, dean_jonction, v_jonction, delta_P_jonction, P2_jonction = \
        coefficient_convection_serpentin(fluid, 1, diam_jonction, None, None, long_jonction, 1.0, \
                                         debit, temp_chaud_out, "maintain")
    if v_jonction < v_min:
        v_min = v_jonction

    # Paramètres
    props = proprietes_fluide(fluid, temp_chaud_out)

    # Calcul
    P_sing = pertes_singulieres(debit, obstacles, props['rho'])
    # Ajouter les coudes et les raccords à la charge des jonctions:
    delta_P_jonction += P_sing

    print("\n"+"="*50)
    print(f"Débit {int(debit*3600*1000+0.5)}L/H:")

    print("\nTempératures:")
    print(f"- Calandre, entrée: {temp_froid_in:.1f}°C, sortie {temp_froid_out_target:.1f}°C, reçoit {q_transfert/1000.0:.3f} kW")
    print(f"- Tube, entrée: {temp_chaud_in:.1f}°C, sortie {temp_chaud_out:.1f}°C, donne {q_transfert/1000.0:.3f} kW")
    print(f"- Appoint: {appoint/1000.0:.3f} kW pour amener {debit * 3600:.3f} m³/h de {temp_froid_out_target:.1f}°C à {temp_chaud_in:.1f} °C")
    print(f"- Jonctions: {long_jonction:.3f} m, diamètre {diam_jonction * mmInM:.3f} mm")

    print("\nConditions opératoires:")
    print(f"- Débit tube: {debit/parallele*3600:.3f} m³/h * {parallele:.1f} paralleles ({fluid} à {temp_tube:.1f}°C)")
    print(f"  Vitesse et Reynold tube: {v_tube:.3f} m/sec, R={rey_tube:.1f}")
    if a_plaque is None:
        print(f"- Débit calandre: {debit/parallele*3600:.3f} m³/h * {parallele:.1f} paralleles ({fluid} à {temp_calandre:.1f}°C)")
        print(f"  Vitesse et Reynold calandre: {v_annulaire:.3f} m/sec, R={rey_annulaire:.1f}")
        # Vitesse maximale dans les coudes
        v_max = v_annulaire * (1 + 0.75 * ( (diam_calandre - diam_ext_tube) / diam_serpentin)**0.5)
        print(f"  Vitesse maximale dans les courbes: {v_max:.3f} m/sec")
    print(f"- Longueur échangeur actuelle: {longueur:.3f} m")
    print(f"  Longueur échangeur calculée selon transfert possible: {longueur2:.3f} m")
    print(f"- Chauffe: {longueur_chauffe:.3f} m, diamètre {diam_chauffe * mmInM:.3f} mm")
    print(f"  Vitesse et Reynold Chauffe: {v_chauffe:.3f} m/sec, R={rey_chauffe:.1f}")
    print(f"- Maintien: {longueur_maintien:.3f} m, diamètre {diam_maintien * mmInM:.3f} mm")
    print(f"  Vitesse et Reynold Maintien: {v_maintien:.3f} m/sec, Durée {longueur_maintien/v_maintien:.1f} sec, R={rey_maintien:.1f}")
    print(f"- Jonctions: {long_jonction:.3f} m, diamètre {diam_jonction * mmInM:.3f} mm")
    print(f"  Vitesse et Reynold Jonctions: {v_jonction:.3f} m/sec, R={rey_jonction:.1f}")

    print("\nÉchangeur:")
    print(f"- LMTD: {lmtd:.1f} °C")
    print(f"- Puissance à échanger: {q_transfert/1000.0:.3f} kW, Efficacité requise selon appoint: {100*q_transfert/(q_transfert+appoint):.1f}%")
    U_necessaire = q_transfert/ lmtd / surface_actuelle
    NTU = (q_transfert / lmtd)/(1000*debit*props['cp'])
    print(f"- Puissance nécessaire par m²K: {U_necessaire:.1f} W/m²K, NTU={NTU:.1f} --> Efficacité: {100*NTU/(1+NTU):.1f}%")
    print(f"  Puissance possible à échanger par m²K: {U_possible:.1f} W/m²K, {U_possible/U_necessaire*100.0:.1f}% du nécessaire")
    print(f"- Surface actuelle: {surface_actuelle:.1f} m²")
    print(f"  Surface calculée selon transfert possible: {surface_echange:.1f} m²")

    print("\nChauffe d'appoint:")
    if diam_calandre_chauffe is None:
        print(f"- Bain marie à {temp_chauffe:.1f}°C")
    else:
        print(f"- Circulation depuis une source chaude à {temp_chauffe:.1f}°C: {debit_circulateur*3600:.1f} m³/h")
    print(f"- Augmentation de température: {temp_chaud_in-temp_froid_out_target:.1f} °C")
    print(f"- Puissance à échanger: {appoint/1000.0:.3f} kW")
    U_chauffe_necessaire = q_transfert/ (temp_chauffe-temp_froid_out_target)/2 / (np.pi * diam_chauffe * longueur_chauffe)
    print(f"- Puissance nécessaire par m²K: {U_chauffe_necessaire:.1f} W/m²K")
    print(f"  Puissance possible à échanger par m²K: {U_possible_chauffe:.1f} W/m²K, {U_possible_chauffe/U_chauffe_necessaire*100.0:.1f}% du nécessaire")

    temp_chauffe_moy = (temp_chaud_in+temp_froid_out_target)/2
    longueurs =    [longueur,     longueur_chauffe,longueur_maintien,    longueur,     long_jonction]
    temperatures = [temp_calandre,temp_chauffe_moy,
                    temp_chaud_in,temp_tube,    temp_calandre]
    if a_plaque:
        diam_plaque = ((largeur_plaque*espacement*parallele) / np.pi)**0.5 * 2
    diametres =    [diam_calandre if a_plaque is None else diam_plaque,
                    diam_chauffe,    diam_maintien,
                    diam_int_tube if a_plaque is None else diam_plaque,   diam_jonction]
    print(f"\nPertes de charge:")
    if a_plaque is None:
        print(f"- Côté froid (annulaire) {temp_calandre:.1f}°C, diamètre {diam_calandre*mmInM:.1f} x {diam_ext_tube*1000:.1f}mm, long={longueur:.1f}m: {doc_P(delta_P_annulaire, P2_annulaire)}")
    else:
        print(f"- Côté froid {temp_calandre:.1f}°C: {doc_P(delta_P_annulaire, P2_annulaire)}")
    # print(f"\n3. Recalcul pour {longueur:.1f} vs {longueur2:.1f} m :")
    # A = np.pi * diam_ext_tube * longueur
    # Q_calcul = u_global * A * (lmtd*longueur2/longueur)
    # print(f"- Puissance thermique potentielle: {Q_calcul/1000:.3f} kW")
    #
    # eau_target = proprietes_fluide(fluid, temp_chaud_in)
    # m_eau = debit_eau * eau_target['rho']  # kg/s
    # cp_eau = eau_target['cp']
    # temp_froid_out_target = temp_froid_in + (Q_calcul / (cp_eau * m_eau))
    # temp_chaud_out  = temp_chaud_in - (temp_froid_out_target - temp_froid_in )
    # print(f"- Calandre, entrée: {temp_froid_in:.1f} °C, sortie {temp_froid_out_target:.1f}°C, reçoit {Q_calcul/1000.0:.3f} kW")
    # print(f"- Tube, entrée: {temp_chaud_in:.1f} °C, sortie {temp_chaud_out:.1f}°C, donne {Q_calcul/1000.0:.3f} kW")
    print(f"- Chauffe {temp_chauffe_moy:.1f}°C, diamètre={diam_chauffe*mmInM:.1f}mm, long={longueur_chauffe:.1f}m: {doc_P(delta_P_chauffe, P2_chauffe)}")

    print(f"- Maintien {duree_maintien:.1f}s  {temp_chaud_in:.1f}°C, diamètre={diam_maintien*mmInM:.1f}mm, long={longueur_maintien:.1f}m: {doc_P(delta_P_maintien,P2_maintien)}")
    if a_plaque is None:
        print(f"- Côté chaud {temp_tube:.1f}°C, diamètre {diam_int_tube*mmInM:.1f}mm, long={longueur:.1f}m: {doc_P(delta_P_tube,P2_tube)}")
    else:
        print(f"- Côté chaud {temp_tube:.1f}°C: {doc_P(delta_P_tube,P2_tube)}")
    print(f"- Jonctions  {temp_calandre:.1f}°C, diamètre={diam_jonction*mmInM:.1f}mm, long={long_jonction:.1f}m: {doc_P(delta_P_jonction,P2_jonction)}")

    total_load = charge_ligne_droite(fluid,debit,longueurs,temperatures,diametres)
    print(f"- TOTAL, long= {np.sum(longueurs):.3f}m:",
            f" {doc_P(delta_P_annulaire+delta_P_tube+delta_P_chauffe+delta_P_maintien+delta_P_jonction,None)}")
    print(f"- en ligne droite  {doc_P(total_load,None)}")
    partJonction = delta_P_jonction/3.0
    results.append([int(debit*1000*3600),(delta_P_annulaire+partJonction)/PaInBar,(delta_P_chauffe+delta_P_maintien+partJonction)/PaInBar,(delta_P_tube+partJonction)/PaInBar])
    #x=f"{(delta_P_annulaire+partJonction)/PaInBar:.3f}\t\t{(delta_P_chauffe+delta_P_maintien+partJonction)/PaInBar:.3f}\t\t{(delta_P_tube+partJonction)/PaInBar:.3f}"
    #print (x.replace('.',','))
    #P2_annulaire+P2_tube+P2_chauffe+P2_maintien)}")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

print(f"\nNETTOYAGE à 1.5 mètres par seconde (Eau)")

if debit_nettoyage is None:
    debit_nettoyage = debit * (v_nettoyage / v_min) # m /s : débit à atteindre pour avoir partout une vitesse de nettoyage correcte
elif v_nettoyage is None:
    v_nettoyage = v_min * (debit_nettoyage / debit)

print (f"- Débit à {v_min:.3f}m/s: {debit*3600000:.1f} L/h, Débit nécessaire pour {v_nettoyage:.3f}m/s: {debit_nettoyage*3600000:.1f} L/h")

if a_plaque is None:
    # """Calcule le coefficient global de transfert thermique U (W/m²K)"""
    # Coefficients de convectionh, f_curved, Re, De, v
    h_tube, f_tube, rey_tube, dean_tube, v_tube, delta_P_tube, P2_tube = \
        coefficient_convection_serpentin(cleaning_fluid, parallele,
                                         diam_int_tube, diam_ext_tube, diam_calandre, longueur, diam_serpentin, \
                                         debit_nettoyage, temp_nettoyage, "tube")
    print(f"- Chaud-->Froid: {doc_P(delta_P_tube,P2_tube)}")
    h_annulaire, f_annulaire, rey_annulaire, dean_annulaire, v_annulaire, delta_P_annulaire, P2_annulaire = \
        coefficient_convection_serpentin(cleaning_fluid, parallele,
                                         diam_int_tube, diam_ext_tube, diam_calandre, longueur, diam_serpentin, \
                                         debit_nettoyage, temp_nettoyage, "shell")
    print(f"- Froid-->Chaud: {doc_P(delta_P_annulaire,P2_annulaire)}")
else:
    h_tube, f_tube, rey_tube, v_tube, delta_P_tube = \
        coefficient_convection_plaque(cleaning_fluid, section_canal, longueur, angle_chevron,
                                      debit_nettoyage, temp_nettoyage, parallele,diam_jonction)
    print(f"- Froid-->Chaud et Chaud-->Froid: {doc_P(delta_P_tube,P2_tube)}")
    delta_P_annulaire = delta_P_tube

h_chauffe, f_chauffe, rey_chauffe, dean_chauffe, v_chauffe, delta_P_chauffe, P2_chauffe = \
    coefficient_convection_serpentin(cleaning_fluid, 1, diam_chauffe, diam_chauffe_ext,
                                     diam_calandre_chauffe, longueur_chauffe, diam_serpentin, \
                                     debit_nettoyage, temp_nettoyage, "tube")
print(f"- Chauffe: {doc_P(delta_P_chauffe, P2_chauffe)}")
h_maintien, f_maintien, rey_maintien, dean_maintien, v_maintien, delta_P_maintien, P2_maintien = \
    coefficient_convection_serpentin(cleaning_fluid, 1, diam_maintien, None, None, longueur_maintien, diam_serpentin, \
                                     debit_nettoyage, temp_nettoyage, "maintain")
print(f"- Maintien: {doc_P(delta_P_maintien,P2_maintien)}")
h_jonction, f_jonction, rey_jonction, dean_jonction, v_jonction, delta_P_jonction, P2_jonction = \
    coefficient_convection_serpentin(fluid, 1, diam_jonction, None, None, long_jonction, 1.0, \
                                     debit, temp_chaud_out, "maintain")
print(f"- Jonction: {doc_P(delta_P_jonction,P2_jonction)}")
print(f"- TOTAL nettoyage, long= {longueur*2+longueur_chauffe+longueur_maintien+long_jonction:.3f}m:",
      f" {doc_P(delta_P_annulaire+delta_P_tube+delta_P_chauffe+delta_P_maintien+delta_P_jonction,None)}")
#P2_annulaire+P2_tube+P2_chauffe+P2_maintien)}")

for i in range(0,len(results)):
    res = results[i]
    print(f"{int(res[0])}\t{res[1]:3f}\t{res[2]:3f}\t{res[3]:3f}")
