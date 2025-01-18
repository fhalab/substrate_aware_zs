"""
This file contains variables specifically to each dataset
"""

# from __future__ import annotations

from copy import deepcopy

# from itertools import combinations

# import numpy as np
# import pandas as pd

# from REVIVAL.util import get_file_name


TrpB_CHEM = {
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
    "substrate-addH_chai_joint": [
        ("B", "LIG", "N1"),
    ],
}

LIB_INFO_DICT = deepcopy(
    {
        "5DW0": {
            "enzyme": "PfTrpB",
            "substrate": "",
            "substrate-smiles": "",
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
        },
        "PfTrpB-4bromo": {
            "enzyme": "PfTrpB",
            "substrate": "4bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C(=C1)Br",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5bromo": {
            "enzyme": "PfTrpB",
            "substrate": "5bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Br",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7bromo": {
            "enzyme": "PfTrpB",
            "substrate": "7bromo",
            "substrate-smiles": "C1=CC2=C(C(=C1)Br)NC=C2",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5chloro": {
            "enzyme": "PfTrpB",
            "substrate": "5chloro",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Cl",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-6chloro": {
            "enzyme": "PfTrpB",
            "substrate": "6chloro",
            "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Cl",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-56chloro": {
            "enzyme": "PfTrpB",
            "substrate": "56chloro",
            "substrate-smiles": "ClC(C=C1NC=CC1=C2)=C2Cl",
            **TrpB_CHEM,
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-4cyano": {
            "enzyme": "PfTrpB",
            "substrate": "4cyano",
            "substrate-smiles": "C1=CC(=C2C=CNC2=C1)C#N",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5cyano": {
            "enzyme": "PfTrpB",
            "substrate": "5cyano",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1C#N",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5iodo": {
            "enzyme": "PfTrpB",
            "substrate": "5iodo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1I",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7iodo": {
            "enzyme": "PfTrpB",
            "substrate": "7iodo",
            "substrate-smiles": "C1=CC2=C(C(=C1)I)NC=C2",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7methyl": {
            "enzyme": "PfTrpB",
            "substrate": "7methyl",
            "substrate-smiles": "CC1=CC=CC2=C1NC=C2",
            **TrpB_CHEM,
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "Rma-CB": {
            "enzyme": "Rma",
            "substrate": "NHC-borane",
            "substrate-smiles": "[BH3-]C1=[N+](C)C=CN1C",
            "cofactor": ["carbene-heme"],
            "cofactor-smiles": [
                "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=C(C)C(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
            ],
            "inactivated-cofactor": ["diazo ester (Me-EDA)", "heme c"],
            "inactivated-cofactor-smiles": [
                "CC(C(OCC)=O)=[N+]=[N-]",
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C",
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
            "substrate-addH_chai_joint": [
                ("B", "LIG", "B1"),
            ],
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
        "Rma-CSi": {
            "enzyme": "Rma",
            "substrate": "phenyldimethylsilane",
            "substrate-smiles": "[H][Si](C)(C)C1=CC=CC=C1",
            "cofactor": ["carbene-heme"],
            "cofactor-smiles": [
                "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=C(C)C(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
            ],
            "inactived-cofactor": ["diazo ester (Me-EDA)", "heme c"],
            "inactivated-cofactor-smiles": [
                "CC(C(OCC)=O)=[N+]=[N-]",
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C",
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
        "ParLQ": {
            "enzyme": "ParLQ",
            "substrate": "4-vinylanisole",
            "substrate-smiles": "C=CC1=CC=C(OC)C=C1",
            "cofactor": ["carbene-heme"],
            "cofactor-smiles": [
                "C=CC1=C(C=C2C(C)=C(C=C)C3=[N]2[Fe]45(N6C(C(C)=C(CCC([O-])=O)C6=C7)=C3)=CC(OCC)=O)N4C(C=C8[N]5=C7C(CCC([O-])=O)=C8C)=C1C",
            ],
            "inactivated-cofactor": ["ethyl diazoacetate (EDA)", "heme b"],
            "inactivated-cofactor-smiles": [
                "[N-]=[N+]=CC(OCC)=O",
                "Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O",
            ],  # heme b taken from pdb
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
            },  # if need to add H
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
        },
    }
)

TmTrpB_LIBS = {
    "TrpB3A": {
        "positions": {1: 104, 2: 105, 3: 106},
        "AAs": {1: "A", 2: "E", 3: "T"},
    },
    "TrpB3B": {
        "positions": {1: 105, 2: 106, 3: 107},
        "AAs": {1: "E", 2: "T", 3: "G"},
    },
    "TrpB3C": {
        "positions": {1: 106, 2: 107, 3: 108},
        "AAs": {1: "T", 2: "G", 3: "A"},
    },
    "TrpB3D": {
        "positions": {1: 117, 2: 118, 3: 119},
        "AAs": {1: "T", 2: "A", 3: "A"},
    },
    "TrpB3E": {
        "positions": {1: 184, 2: 185, 3: 186},
        "AAs": {1: "F", 2: "G", 3: "S"},
    },
    "TrpB3F": {
        "positions": {1: 162, 2: 166, 3: 301},
        "AAs": {1: "L", 2: "I", 3: "Y"},
    },
    "TrpB3G": {
        "positions": {1: 227, 2: 228, 3: 301},
        "AAs": {1: "V", 2: "S", 3: "Y"},
    },
    "TrpB3H": {
        "positions": {1: 228, 2: 230, 3: 231},
        "AAs": {1: "S", 2: "G", 3: "S"},
    },
    "TrpB3I": {
        "positions": {1: 182, 2: 183, 3: 184},
        "AAs": {1: "Y", 2: "V", 3: "F"},
    },
    "TrpB4": {
        "positions": {1: 183, 2: 184, 3: 227, 4: 228},
        "AAs": {1: "V", 2: "F", 3: "V", 4: "S"},
    },
}

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