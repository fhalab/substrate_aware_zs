"""
This file contains variables specifically to each dataset
"""

from copy import deepcopy


AA_DICT = {
    "A": "ALA",  # Alanine
    "R": "ARG",  # Arginine
    "N": "ASN",  # Asparagine
    "D": "ASP",  # Aspartic acid
    "C": "CYS",  # Cysteine
    "E": "GLU",  # Glutamic acid
    "Q": "GLN",  # Glutamine
    "G": "GLY",  # Glycine
    "H": "HIS",  # Histidine
    "I": "ILE",  # Isoleucine
    "L": "LEU",  # Leucine
    "K": "LYS",  # Lysine
    "M": "MET",  # Methionine
    "F": "PHE",  # Phenylalanine
    "P": "PRO",  # Proline
    "S": "SER",  # Serine
    "T": "THR",  # Threonine
    "W": "TRP",  # Tryptophan
    "Y": "TYR",  # Tyrosine
    "V": "VAL",  # Valine
}

TrpB_COMMON = {
    "cofactor": ["PLP-dependent_aminoacrylate", "Na+"],
    "cofactor-smiles": [
        "[O-]C1=C(/C=[N+]([H])/C(C([O-])=O)=C)C(CP([O-])([O-])=O)=CN=C1C",
        "[Na+]",
    ],
    "inactivated-cofactor": ["PLP", "Na+"],
    "inactivated-cofactor-smiles": [
        "O=Cc1c(O)c(C)ncc1COP(O)(O)=O",
        "[Na+]",
    ],
    "transiminated-cofactor": ["aldimine", "Na+"],
    "transiminated-cofactor-smiles": [
        "[O-]C1=C(/C=[N+]([H])/[C@](CO)([H])C([O-])=O)C(CP([O-])([O-])=O)=CN=C1C",
        "[Na+]",
    ],
    "positions": {1: 165, 2: 183, 3: 301},
    "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
    "family": "TrpB",
    "project": "multi-substrate",
}

ParLQ_COMMON = {
    "cofactor": ["carbene-heme"],
    "cofactor-smiles": [
        "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=CC(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
    ],
    "carbene_precursor": "ethyl diazoacetate (EDA)",
    "carbene_precursor-smiles": "[N-]=[N+]=CC(OCC)=O",
    "inactivated-cofactor": ["heme b"],
    "inactivated-cofactor-smiles": [
        "Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O",
    ],  # heme b taken from pdb
    "activated_carbene-cofactor": ["heme c", "Fe2+"," carbene"],
    "activated_carbene-cofactor-smiles": [
        r"C=CC1=C(/C=C2C(C)=C(C=C)C3=N/2)NC(/C=C4N=C(/C=C(C(CCC([O-])=O)=C/5C)\NC5=C/3)C(CCC([O-])=O)=C\4C)=C1C",
        "[Fe2+]",
        "CCOC([C])=O"
    ],
    "cofactor-distances_seperate": {
        "C-C_1": (
            ("C", 1, "LIG", "C_15", False),
            ("B", 1, "LIG", "C_1", False),
        ),
        "C-C_2": (
            ("C", 1, "LIG", "C_15", False),
            ("B", 1, "LIG", "C_2", False),
        ),
    },  # if need to add H
    "positions": {1: 56, 2: 57, 3: 59, 4: 60, 5: 89},
    "AAs": {
        1: "W",
        2: "Y",
        3: "L",
        4: "Q",
        5: "F",
    },  # W56, Y57, L59, Q60, and F89; WYLQF
    "family": "ParPgb",
    "project": "ALDE",
}

ENZYME_INFO_DICT = deepcopy(
    {
        "PfTrpB": {
            "ligand-info": [
                ("A", "PLS", "N"),
                ("A", "PLS", "CA"),
                ("A", "PLS", "CB"),
                ("A", "PLS", "OG"),
                ("A", "PLS", "C"),
                ("A", "PLS", "O"),
                ("A", "PLS", "OXT"),
                ("A", "PLS", "N1"),
                ("A", "PLS", "C2"),
                ("A", "PLS", "C2A"),
                ("A", "PLS", "C3"),
                ("A", "PLS", "O3"),
                ("A", "PLS", "C4"),
                ("A", "PLS", "C4A"),
                ("A", "PLS", "C5"),
                ("A", "PLS", "C6"),
                ("A", "PLS", "C5A"),
                ("A", "PLS", "O4P"),
                ("A", "PLS", "P"),
                ("A", "PLS", "O1P"),
                ("A", "PLS", "O2P"),
                ("A", "PLS", "O3P"),
            ],
            "inactivated-cofactor": ["PLP", "Na+"],
            "inactivated-cofactor-smiles": [
                "O=Cc1c(O)c(C)ncc1COP(O)(O)=O",
                "[Na+]",
            ],
            "transiminated-cofactor": ["aldimine", "Na+"],
            "transiminated-cofactor-smiles": [
                "[O-]C1=C(/C=[N+]([H])/[C@](CO)([H])C([O-])=O)C(CP([O-])([O-])=O)=CN=C1C",
                "[Na+]",
            ],
            "cofactor": ["PLS", "Na+"],
            "cofactor-smiles": ["CC1=NC=C(COP(O)(O)=O)C(CNC(C(O)=O)CO)=C1O", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
            "PDBID": "5DW0",
            "volume": 792.872,  # CASTp with 5dw3 PocID2
        },
        "Rma": {
            "ligand-info": [
                ("A", "SO4", "S"),
                ("A", "SO4", "O1"),
                ("A", "SO4", "O2"),
                ("A", "SO4", "O3"),
                ("A", "SO4", "O4"),
                ("A", "HEC", "FE"),
                ("A", "HEC", "CHA"),
                ("A", "HEC", "CHB"),
                ("A", "HEC", "CHC"),
                ("A", "HEC", "CHD"),
                ("A", "HEC", "NA"),
                ("A", "HEC", "C1A"),
                ("A", "HEC", "C2A"),
                ("A", "HEC", "C3A"),
                ("A", "HEC", "C4A"),
                ("A", "HEC", "CMA"),
                ("A", "HEC", "CAA"),
                ("A", "HEC", "CBA"),
                ("A", "HEC", "CGA"),
                ("A", "HEC", "O1A"),
                ("A", "HEC", "O2A"),
                ("A", "HEC", "NB"),
                ("A", "HEC", "C1B"),
                ("A", "HEC", "C2B"),
                ("A", "HEC", "C3B"),
                ("A", "HEC", "C4B"),
                ("A", "HEC", "CMB"),
                ("A", "HEC", "CAB"),
                ("A", "HEC", "CBB"),
                ("A", "HEC", "NC"),
                ("A", "HEC", "C1C"),
                ("A", "HEC", "C2C"),
                ("A", "HEC", "C3C"),
                ("A", "HEC", "C4C"),
                ("A", "HEC", "CMC"),
                ("A", "HEC", "CAC"),
                ("A", "HEC", "CBC"),
                ("A", "HEC", "ND"),
                ("A", "HEC", "C1D"),
                ("A", "HEC", "C2D"),
                ("A", "HEC", "C3D"),
                ("A", "HEC", "C4D"),
                ("A", "HEC", "CMD"),
                ("A", "HEC", "CAD"),
                ("A", "HEC", "CBD"),
                ("A", "HEC", "CGD"),
                ("A", "HEC", "O1D"),
                ("A", "HEC", "O2D"),
            ],
            "family": "Rma",
            "project": "MODIFY",
            "PDBID": "3CP5",
            "volume": 138.961,  # CASTp
        },
        "ParLQ": {
            "ligand-info": [
                ("B", "HEM", "CAA"),
                ("B", "HEM", "CAB"),
                ("B", "HEM", "CAC"),
                ("B", "HEM", "CAD"),
                ("B", "HEM", "NA"),
                ("B", "HEM", "CBA"),
                ("B", "HEM", "CBB"),
                ("B", "HEM", "CBC"),
                ("B", "HEM", "CBD"),
                ("B", "HEM", "NB"),
                ("B", "HEM", "CGA"),
                ("B", "HEM", "CGD"),
                ("B", "HEM", "ND"),
                ("B", "HEM", "CHA"),
                ("B", "HEM", "CHB"),
                ("B", "HEM", "CHC"),
                ("B", "HEM", "CHD"),
                ("B", "HEM", "CMA"),
                ("B", "HEM", "CMB"),
                ("B", "HEM", "CMC"),
                ("B", "HEM", "CMD"),
                ("B", "HEM", "C1A"),
                ("B", "HEM", "C1B"),
                ("B", "HEM", "C1C"),
                ("B", "HEM", "C1D"),
                ("B", "HEM", "O1A"),
                ("B", "HEM", "O1D"),
                ("B", "HEM", "C2A"),
                ("B", "HEM", "C2B"),
                ("B", "HEM", "C2C"),
                ("B", "HEM", "C2D"),
                ("B", "HEM", "O2A"),
                ("B", "HEM", "O2D"),
                ("B", "HEM", "C3A"),
                ("B", "HEM", "C3B"),
                ("B", "HEM", "C3C"),
                ("B", "HEM", "C3D"),
                ("B", "HEM", "C4A"),
                ("B", "HEM", "C4B"),
                ("B", "HEM", "C4C"),
                ("B", "HEM", "C4D"),
                ("B", "HEM", "NAC"),
                ("B", "HEM", "FE"),
            ],
            "positions": {1: 56, 2: 57, 3: 59, 4: 60, 5: 89},
            "AAs": {
                1: "W",
                2: "Y",
                3: "L",
                4: "Q",
                5: "F",
            },  # W56, Y57, L59, Q60, and F89; WYLQF\
            "family": "ParLQ",
            "project": "ALDE",
            "volume": 275.468,  # CASTp with 3zji
        },
    }
)

LIB_INFO_DICT = deepcopy(
    {
        "PfTrpB-4bromo": {
            "enzyme": "PfTrpB",
            "substrate": "4bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C(=C1)Br",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(Br)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "4bromo-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "BR1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C7", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-5bromo": {
            "enzyme": "PfTrpB",
            "substrate": "5bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Br",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Br)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "5bromo-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "BR1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C6", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
        },
        "PfTrpB-7bromo": {
            "enzyme": "PfTrpB",
            "substrate": "7bromo",
            "substrate-smiles": "C1=CC2=C(C(=C1)Br)NC=C2",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(Br)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "7bromo-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "BR1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_8", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C6", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-5chloro": {
            "enzyme": "PfTrpB",
            "substrate": "5chloro",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Cl",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "5chloro-info": [
                ("B", "LIG", "CL1"),
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C6", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-6chloro": {
            "enzyme": "PfTrpB",
            "substrate": "6chloro",
            "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Cl",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "6chloro-info": [
                ("B", "LIG", "CL1"),
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_7", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C5", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-56chloro": {
            "enzyme": "PfTrpB",
            "substrate": "56chloro",
            "substrate-smiles": "ClC(C=C1NC=CC1=C2)=C2Cl",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "56chloro-info": [
                ("B", "LIG", "CL1"),
                ("B", "LIG", "CL2"),
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C4", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-4cyano": {
            "enzyme": "PfTrpB",
            "substrate": "4cyano",
            "substrate-smiles": "C1=CC(=C2C=CNC2=C1)C#N",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(C#N)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "4cyano-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "C9"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "N2"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N3"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N3"), ("B", "LIG", "C12")),
                (("B", "LIG", "C13"), ("B", "LIG", "C15")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_15", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N2"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N3"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N3"), ("B", "LIG", "C12")),
                (("B", "LIG", "C10"), ("B", "LIG", "C11")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C8", False), ("B", 1, "LIG", "C10", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N2", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N2", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-5cyano": {
            "enzyme": "PfTrpB",
            "substrate": "5cyano",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1C#N",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(C=C23)C#N)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "5cyano-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "C9"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "N2"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N3"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N3"), ("B", "LIG", "C12")),
                (("B", "LIG", "C13"), ("B", "LIG", "C15")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_15", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N2"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N3"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N3"), ("B", "LIG", "C12")),
                (("B", "LIG", "C10"), ("B", "LIG", "C11")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C7", False), ("B", 1, "LIG", "C10", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N2", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N2", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-5iodo": {
            "enzyme": "PfTrpB",
            "substrate": "5iodo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1I",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(I)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "5iodo-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "I1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_5", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C6", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
        #     "positions": {1: 165, 2: 183, 3: 301},
        #     "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
        #     "family": "TrpB",
        #     "project": "multi-substrate",
        },
        "PfTrpB-7iodo": {
            "enzyme": "PfTrpB",
            "substrate": "7iodo",
            "substrate-smiles": "C1=CC2=C(C(=C1)I)NC=C2",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(I)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "7iodo-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "I1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C12"), ("B", "LIG", "C14")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_8", False), ("B", 1, "LIG", "C_14", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O6"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C11")),
                (("B", "LIG", "C9"), ("B", "LIG", "C10")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C6", False), ("B", 1, "LIG", "C9", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB-7methyl": {
            "enzyme": "PfTrpB",
            "substrate": "7methyl",
            "substrate-smiles": "CC1=CC=CC2=C1NC=C2",
            **TrpB_COMMON,
            "intermediate-smiles": [
                "CC1=C2NC=C(C[C@H](\[NH+]=C\C3=C(COP([O-])([O-])=O)C=NC(C)=C3[O-])C([O-])=O)C2=CC=C1",
                "[Na+]",
            ],
            "7methyl-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "C9"),
                ("B", "LIG", "N1"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_chai_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O4"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_chai_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C12")),
                (("B", "LIG", "C13"), ("B", "LIG", "C15")),
            ],
            "cofactor-distances_chai_joint": {
                "C-C": (("B", 1, "LIG", "C_9", False), ("B", 1, "LIG", "C_15", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N_1", True),
                ),
            },  # if need to add H
            "substrate-addH_af3": [
                ("B", "LIG", "N1"),
            ],
            "cofactor-addH_af3_joint": [
                ("B", "LIG", "N2"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
                ("B", "LIG", "O5"),
            ],
            "cofactor-double_bond_pairs_af3_joint": [
                (("B", "LIG", "N2"), ("B", "LIG", "C12")),
                (("B", "LIG", "C10"), ("B", "LIG", "C11")),
            ],
            "cofactor-distances_af3_joint": {
                "C-C": (("B", 1, "LIG", "C7", False), ("B", 1, "LIG", "C10", False)),
                "GLU-NH_1": (
                    ("A", 104, "GLU", "OE1", False),
                    ("B", 1, "LIG", "N1", True),
                ),
                "GLU-NH_2": (
                    ("A", 104, "GLU", "OE2", False),
                    ("B", 1, "LIG", "N1", True),
                ),
            },  # if need to add H
            # "positions": {1: 165, 2: 183, 3: 301},
            # "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            # "family": "TrpB",
            # "project": "multi-substrate",
        },
        "PfTrpB_scope":{
            "enzyme": "PfTrpB",
            "positions":{1: 104, 2: 139, 3: 161, 4: 165, 5: 166, 6: 183, 7: 186, 8: 212, 9: 301},
            "AAs":{1: "E", 2: "L", 3: "L", 4: "I", 5: "D", 6: "I", 7: "V", 8: "P", 9: "Y"},
            "family": "TrpB",
            "project": "substrate-scope",
        },
        "Rma-CB": {
            "enzyme": "Rma",
            "substrate": "NHC-borane",
            "substrate-smiles": "[BH3-]C1=[N+](C)C=CN1C",
            "cofactor": ["carbene-heme"],
            "cofactor-smiles": [
                "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=C(C)C(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
            ],
            "carbene_precursor": "diazo ester (Me-EDA)",
            "carbene_precursor-smiles": "CC(C(OCC)=O)=[N+]=[N-]",
            "inactivated-cofactor": ["heme c"],
            "inactivated-cofactor-smiles": [
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C"
            ],
            "activated_carbene-cofactor": ["heme c", "Fe2+", "carbene"],
            "activated_carbene-cofactor-smiles": [
                r"C=CC1=C(/C=C2C(C)=C(C=C)C3=N/2)NC(/C=C4N=C(/C=C(C(CCC([O-])=O)=C/5C)\NC5=C/3)C(CCC([O-])=O)=C\4C)=C1C",
                "[Fe2+]",
                "CCOC([C]C)=O",
            ],
            "product": "organoborane",
            "product-smiles": "CN1C=C[N+](C)=C1[BH2-][C@H](C)C(OCC)=O",
            "NHC-borane-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "N1"),
                ("B", "LIG", "N2"),
                ("B", "LIG", "B1"),
            ],
            "carbene-info_joint": [
                ("B", "LIG", "C20"),
                ("B", "LIG", "C21"),
                ("B", "LIG", "C22"),
                ("B", "LIG", "C23"),
                ("B", "LIG", "C24"),
                ("B", "LIG", "O1"),
                ("B", "LIG", "O2"),
            ],
            "carbene-info_seperate": [
                ("C", "LIG", "C15"),
                ("C", "LIG", "C16"),
                ("C", "LIG", "C17"),
                ("C", "LIG", "C18"),
                ("C", "LIG", "C19"),
                ("C", "LIG", "O1"),
                ("C", "LIG", "O2"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "B1"),
            ],
            "substrate-double_bond_pairs_chai": [
                (("B", "LIG", "C1"), ("B", "LIG", "N1")),
                (("B", "LIG", "C3"), ("B", "LIG", "C4")),
            ],
            "substrate_branches_chai": {"B": "B1", "C": "C1"},
            "substrate-addH_af3": [
                ("B", "LIG", "B1"),
            ],
            "substrate-double_bond_pairs_af3": [
                (("B", "LIG", "C1"), ("B", "LIG", "N1")),
                (("B", "LIG", "C3"), ("B", "LIG", "C4")),
            ],
            "substrate_branches_af3": {"B": "B1", "C": "C1"},
            "cofactor-distances_joint": {
                "C-B": (("B", 1, "LIG", "C_20", False), ("B", 1, "LIG", "B_1", False))
            },  # if need to add H
            "cofactor-distances_seperate": {
                "C-B": (("C", 1, "LIG", "C_15", False), ("B", 1, "LIG", "B_1", False))
            },  # if need to add H
            "positions": {1: 75, 2: 99, 3: 100, 4: 101, 5: 102, 6: 103},
            "AAs": {1: "V", 2: "M", 3: "M", 4: "T", 5: "D", 6: "M"},
            "family": "cyt c",
            "project": "MODIFY",
        },
        "Rma-CB_scope":{
            "enzyme": "Rma",
            "positions": {1: 71, 2: 75, 3: 89, 4: 98, 5: 99, 6: 100, 7: 101, 8: 103},
            "AAs": {1: "Y", 2: "V", 3: "M", 4: "T", 5: "M", 6: "M", 7: "T", 8: "M"},
            "family": "cyt c",
            "project": "substrate-scope",
        },
        "Rma-CSi": {
            "enzyme": "Rma",
            "substrate": "phenyldimethylsilane",
            "substrate-smiles": "[H][Si](C)(C)C1=CC=CC=C1",
            "cofactor": ["carbene-heme"],
            "cofactor-smiles": [
                "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=C(C)C(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
            ],
            "carbene_precursor": "diazo ester (Me-EDA)",
            "carbene_precursor-smiles": "CC(C(OCC)=O)=[N+]=[N-]",
            "inactivated-cofactor": ["heme c"],
            "inactivated-cofactor-smiles": [
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C"
            ],
            "activated_carbene-cofactor": ["heme c", "Fe2+", "carbene"],
            "activated_carbene-cofactor-smiles": [
                r"C=CC1=C(/C=C2C(C)=C(C=C)C3=N/2)NC(/C=C4N=C(/C=C(C(CCC([O-])=O)=C/5C)\NC5=C/3)C(CCC([O-])=O)=C\4C)=C1C",
                "[Fe2+]",
                "CCOC([C]C)=O"
            ],
            "product": "organosilicon",
            "product-smiles": "C[Si]([C@]([H])(C)C(OCC)=O)(C)C1=CC=CC=C1",
            "phenyldimethylsilane-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "SI1"),
            ],
            "carbene-info_joint": [
                ("B", "LIG", "C23"),
                ("B", "LIG", "C24"),
                ("B", "LIG", "C25"),
                ("B", "LIG", "C26"),
                ("B", "LIG", "C27"),
                ("B", "LIG", "O1"),
                ("B", "LIG", "O2"),
            ],
            "carbene-info_seperate": [
                ("C", "LIG", "C15"),
                ("C", "LIG", "C16"),
                ("C", "LIG", "C17"),
                ("C", "LIG", "C18"),
                ("C", "LIG", "C19"),
                ("C", "LIG", "O1"),
                ("C", "LIG", "O2"),
            ],
            "substrate-addH_chai": [
                ("B", "LIG", "SI1"),
            ],
            "substrate-addH_af3": [
                ("B", "LIG", "SI1"),
            ],
            "cofactor-distances_joint": {
                "C-Si": (("B", 1, "LIG", "C_23", False), ("B", 1, "LIG", "SI_1", False))
            },  # if need to add H
            "cofactor-distances_seperate": {
                "C-Si": (("C", 1, "LIG", "C_15", False), ("B", 1, "LIG", "SI_1", False))
            },  # if need to add H
            "positions": {1: 75, 2: 99, 3: 100, 4: 101, 5: 102, 6: 103},
            "AAs": {1: "V", 2: "M", 3: "M", 4: "T", 5: "D", 6: "M"},
            "family": "cyt c",
            "project": "MODIFY",
        },
        "Rma-CSi_scope":{
            "enzyme": "Rma",
            "positions": {1: 75, 2: 100, 3: 103},
            "AAs": {1: "V", 2: "M", 3: "M"},
            "family": "cyt c",
            "project": "substrate-scope",
        },
        "ParLQ": {
            "enzyme": "ParLQ",
            "substrate": "4-vinylanisole",
            "substrate-smiles": "C=CC1=CC=C(OC)C=C1",
            "product": "1,2-disubstituted cyclopropanes cis",
            "product-smiles": "COC1=CC=C([C@@H]2[C@H](C(OCC)=O)C2)C=C1",
            "product2": "1,2-disubstituted cyclopropanes trans",
            "product2-smiles": "COC1=CC=C([C@@H]2[C@@H](C(OCC)=O)C2)C=C1",
            "4-vinylanisole-info": [
                ("B", "LIG", "C1"),
                ("B", "LIG", "C2"),
                ("B", "LIG", "C3"),
                ("B", "LIG", "C4"),
                ("B", "LIG", "C5"),
                ("B", "LIG", "C6"),
                ("B", "LIG", "C7"),
                ("B", "LIG", "C8"),
                ("B", "LIG", "C9"),
                ("B", "LIG", "O1"),
            ],
            "carbene-info_joint": [
                ("B", "LIG", "C24"),
                ("B", "LIG", "C25"),
                ("B", "LIG", "C26"),
                ("B", "LIG", "C27"),
                ("B", "LIG", "O2"),
                ("B", "LIG", "O3"),
            ],
            "carbene-info_seperate": [
                ("C", "LIG", "C15"),
                ("C", "LIG", "C16"),
                ("C", "LIG", "C17"),
                ("C", "LIG", "C18"),
                ("C", "LIG", "O1"),
                ("C", "LIG", "O2"),
            ],
            "cofactor-distances_joint": {
                "C-C_1": (
                    ("B", 1, "LIG", "C_24", False),
                    ("B", 1, "LIG", "C_1", False),
                ),
                "C-C_2": (
                    ("B", 1, "LIG", "C_24", False),
                    ("B", 1, "LIG", "C_2", False),
                ),
            }, 
            **ParLQ_COMMON,
        },
        "ParLQ-b": {
            "enzyme": "ParLQ",
            "substrate": "b",
            "substrate-smiles": "C=CC1=CC=CC=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-c": {
            "enzyme": "ParLQ",
            "substrate": "c",
            "substrate-smiles": "CC1=CC=C(C=C)C=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-d": {
            "enzyme": "ParLQ",
            "substrate": "d",
            "substrate-smiles": "C=CC1=CC(C)=CC=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-e": {
            "enzyme": "ParLQ",
            "substrate": "e",
            "substrate-smiles": "C=CC1=C(C)C=CC=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-f": {
            "enzyme": "ParLQ",
            "substrate": "f",
            "substrate-smiles": "ClC1=CC=C(C=C)C=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-g": {
            "enzyme": "ParLQ",
            "substrate": "g",
            "substrate-smiles": "BrC1=CC=C(C=C)C=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-h": {
            "enzyme": "ParLQ",
            "substrate": "h",
            "substrate-smiles": "C=CC1=CC=C(C(F)(F)F)C=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
        "ParLQ-i": {
            "enzyme": "ParLQ",
            "substrate": "i",
            "substrate-smiles": "C=CC1=CC(C=CC=C2)=C2C=C1",
            "product": "",
            "product-smiles": "",
            "product2": "",
            "product2-smiles": "",
            **ParLQ_COMMON,
        },
    }
)

APPEND_INFO_COLS = [
    "enzyme",
    "substrate",
    "substrate-smiles",
    "cofactor",
    "cofactor-smiles",
    "intermediate-smiles",
    "product",
    "product-smiles",
]