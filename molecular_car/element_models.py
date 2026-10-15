"""Versioned race support, independent of the molecular and Au force fields.

Membership is an admission boundary, not a claim that every valence or charge
state is parameterized. RDKit UFF and strict Gasteiger checks remain mandatory.
"""

EFFECTIVE_POLARIZABILITY = {
    "H": .20, "B": 1.20, "C": 1.00, "N": .85, "O": .70, "F": .40,
    "Si": 2.50, "P": 1.80, "S": 1.60, "Cl": 1.30, "Br": 1.80, "I": 2.80,
}
SUPPORTED_RACE_ELEMENTS = tuple(EFFECTIVE_POLARIZABILITY)
PROVISIONAL_AU_ELEMENTS = tuple(e for e in SUPPORTED_RACE_ELEMENTS if e not in ("C", "H"))
ELEMENT_MODEL_VERSION = "organic-12-uff-gasteiger-effective-response-v1"


def element_support_metadata(elements, surface_mode="periodic-au"):
    used = sorted(set(elements))
    return {
        "version": ELEMENT_MODEL_VERSION,
        "allowed_elements": list(SUPPORTED_RACE_ELEMENTS),
        "used_elements": used,
        "molecular_parameters": "RDKit UFF, checked for the complete molecular graph",
        "charge_parameters": "strict Gasteiger; unavailable environments rejected",
        "response_parameters": "fixed scalar game coefficients, not measured atomic polarizabilities",
        "effective_response_coefficients": {e: EFFECTIVE_POLARIZABILITY[e] for e in used},
        "gold_interface": {e: "effective element-independent Morse/triangular game surface" if surface_mode == "continuum"
                           else "CHAu SAIP with bonded game normals" if e in ("C", "H")
                           else "provisional nonreactive UFF-vdW mixed lattice LJ" for e in used},
        "limitations": ["Element support does not imply all oxidation states or coordination environments.",
                        "No Au chemisorption, charge transfer, image-charge response or bond breaking."],
    }
