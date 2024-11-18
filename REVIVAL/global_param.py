"""
This file contains variables specifically to each dataset
"""

from __future__ import annotations

from copy import deepcopy
from itertools import combinations

import numpy as np
import pandas as pd

from REVIVAL.util import get_file_name

LIB_INFO_DICT = deepcopy(
    {
        "5DW0": {
            "enzyme": "PfTrpB",
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
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(Br)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5bromo": {
            "enzyme": "PfTrpB",
            "substrate": "5bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Br",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Br)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        # "PfTrpB-6bromo": {
        #     "enzyme": "PfTrpB",
        #     "substrate": "6bromo",
        #     "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Br",
        #     "cofactor": ["PLP", "Na+"],
        #     "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
        #     "positions": {1: 165, 2: 183, 3: 301},
        #     "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
        #     "family": "TrpB",
        #     "project": "multi-substrate",
        # },
        "PfTrpB-7bromo": {
            "enzyme": "PfTrpB",
            "substrate": "7bromo",
            "substrate-smiles": "C1=CC2=C(C(=C1)Br)NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(Br)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5chloro": {
            "enzyme": "PfTrpB",
            "substrate": "5chloro",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-6chloro": {
            "enzyme": "PfTrpB",
            "substrate": "6chloro",
            "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-56chloro": {
            "enzyme": "PfTrpB",
            "substrate": "56chloro",
            "substrate-smiles": "ClC(C=C1NC=CC1=C2)=C2Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-4cyano": {
            "enzyme": "PfTrpB",
            "substrate": "4cyano",
            "substrate-smiles": "C1=CC(=C2C=CNC2=C1)C#N",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(C#N)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5cyano": {
            "enzyme": "PfTrpB",
            "substrate": "5cyano",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1C#N",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(C=C23)C#N)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5iodo": {
            "enzyme": "PfTrpB",
            "substrate": "5iodo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1I",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(I)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7iodo": {
            "enzyme": "PfTrpB",
            "substrate": "7iodo",
            "substrate-smiles": "C1=CC2=C(C(=C1)I)NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(I)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7methyl": {
            "enzyme": "PfTrpB",
            "substrate": "7methyl",
            "substrate-smiles": "CC1=CC=CC2=C1NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "intermediate-smiles": [
                "CC1=C2NC=C(C[C@H](\[NH+]=C\C3=C(COP([O-])([O-])=O)C=NC(C)=C3[O-])C([O-])=O)C2=CC=C1",
                "[Na+]",
            ],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"},  # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "Rma-CB": {
            "enzyme": "Rma",
            "substrate": "NHC-borane",
            "substrate-smiles": "[BH3-]C1=[N+](C)C=CN1C",
            "cofactor": ["diazo ester (Me-EDA)", "heme c"],
            "cofactor-smiles": [
                "CC(C(OCC)=O)=[N+]=[N-]",
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C",
            ],
            "product": "organoborane",
            "product-smiles": "CN1C=C[N+](C)=C1[BH2-][C@H](C)C(OCC)=O",
            "positions": {1: 75, 2: 99, 3: 100, 4: 101, 5: 102, 6: 103},
            "AAs": {1: "V", 2: "M", 3: "M", 4: "T", 5: "D", 6: "M"},
            "family": "cyt c",
            "project": "MODIFY",
        },
        "Rma-CSi": {
            "enzyme": "Rma",
            "substrate": "phenyldimethylsilane",
            "substrate-smiles": "[H][Si](C)(C)C1=CC=CC=C1",
            "cofactor": ["diazo ester (Me-EDA)", "heme c"],
            "cofactor-smiles": [
                "CC(C(OCC)=O)=[N+]=[N-]",
                "CC=C1C(=C2C=C3C(=CC)C(=C4N3[Fe]56N2C1=Cc7n5c(c(c7C)CCC(=O)O)C=C8N6C(=C4)C(=C8CCC(=O)O)C)C)C",
            ],
            "product": "organosilicon",
            "product-smiles": "C[Si]([C@]([H])(C)C(OCC)=O)(C)C1=CC=CC=C1",
            "positions": {1: 75, 2: 99, 3: 100, 4: 101, 5: 102, 6: 103},
            "AAs": {1: "V", 2: "M", 3: "M", 4: "T", 5: "D", 6: "M"},
            "family": "cyt c",
            "project": "MODIFY",
        },
        "ParLQ": {
            "enzyme": "ParLQ",
            "substrate": "4-vinylanisole",
            "substrate-smiles": "C=CC1=CC=C(OC)C=C1",
            "cofactor": ["ethyl diazoacetate (EDA)", "heme b"],
            "cofactor-smiles": [
                "[N-]=[N+]=CC(OCC)=O",
                "Cc1c2n3c(c1CCC(=O)O)C=C4C(=C(C5=[N]4[Fe]36[N]7=C(C=C8N6C(=C5)C(=C8C)C=C)C(=C(C7=C2)C)C=C)C)CCC(=O)O"
            ], # heme b taken from pdb
            "product": "1,2-disubstituted cyclopropanes cis",
            "product-smiles": "COC1=CC=C([C@@H]2[C@H](C(OCC)=O)C2)C=C1",
            "positions": {1: 56, 2: 57, 3: 59, 4: 60, 5: 89},
            "AAs": {1: "W", 2: "Y", 3: "L", 4: "Q", 5: "F"}, # W56, Y57, L59, Q60, and F89; WYLQF
            "family": "ParPgb",
            "project": "ALDE",
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