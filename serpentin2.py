import math
PaInBar = 100000
PaInmBar = 100

def coefficients_tube_en_serpentin(
        fluid: str,
        diam_int_tube: float,
        longueur: float,
        diam_serpentin: float,
        debit: float,
        temp: float
) -> dict:
    """
    Calcule la perte de charge dans un tube intérieur en serpentin (forme hélicoïdale).

    Paramètres :
    - fluid : 'water', 'milk'
    - diam_int_tube : diamètre intérieur du tube (m)
    - longueur : longueur totale développée du serpentin (m)
    - diam_serpentin : diamètre moyen de la spirale du serpentin (m)
    - debit : débit volumique en m³/h
    - temp : température du fluide (°C)

    Retourne :
    - f : facteur de friction corrigé Dean
    - Re : Reynolds
    - De : Dean
    - v : vitesse moyenne (m/s)
    - delta_P : perte de charge totale (Pa)
    """

    # Propriétés du fluide
    if fluid.lower() == "milk":
        rho = 1020  # kg/m³
        mu = 1.9e-3  # Pa.s à 85°C
    elif fluid.lower() == "water":
        rho = 970
        mu = 0.36e-3
    else:
        raise ValueError("Fluide non pris en charge")

    Q_m3_s = debit
    A = math.pi * (diam_int_tube**2) / 4
    v = Q_m3_s / A
    Re = rho * v * diam_int_tube / mu

    # Nombre de Dean
    De = Re * math.sqrt(diam_int_tube / diam_serpentin)

    # Facteur de correction empirique selon Dean (White, 2003)
    # f_corrigé = f_droit * (1 + 0.033 * De)
    if Re < 2300:
        f_droit = 64 / Re
    else:
        f_droit = 0.3164 / (Re**0.25)

    f_corrige = f_droit * min(1 + 0.033 * De, 3.0)  # validé expérimentalement pour De < 500

    delta_P = f_corrige * (longueur / diam_int_tube) * (rho * v**2 / 2)

    return {
        "f": f_corrige,
        "Reynold": Re,
        "Dean": De,
        "vitesse": v,
        "delta_P": delta_P / PaInmBar
    }

# Exemple d'appel
resultats = coefficients_tube_en_serpentin(
    fluid="water",
    diam_int_tube=0.0095,
    longueur=18.0,
    diam_serpentin=0.35,
    debit=0.250/3600,
    temp=85
)

for cle, valeur in resultats.items():
    print(f"{cle} : {valeur:.3f}")
