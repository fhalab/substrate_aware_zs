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
        "PfTrpB-4bromo": {
            "enzyme": "PfTrpB",
            "substrate": "4bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C(=C1)Br",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(Br)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5bromo": {
            "enzyme": "PfTrpB",
            "substrate": "5bromo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Br",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Br)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
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
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(Br)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5chloro": {
            "enzyme": "PfTrpB",
            "substrate": "5chloro",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-6chloro": {
            "enzyme": "PfTrpB",
            "substrate": "6chloro",
            "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-56chloro": {
            "enzyme": "PfTrpB",
            "substrate": "56chloro",
            "substrate-smiles": "ClC(C=C1NC=CC1=C2)=C2Cl",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC(Cl)=C(Cl)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-4cyano": {
            "enzyme": "PfTrpB",
            "substrate": "4cyano",
            "substrate-smiles": "C1=CC(=C2C=CNC2=C1)C#N",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=CC(C#N)=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5cyano": {
            "enzyme": "PfTrpB",
            "substrate": "5cyano",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1C#N",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(C=C23)C#N)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-5iodo": {
            "enzyme": "PfTrpB",
            "substrate": "5iodo",
            "substrate-smiles": "C1=CC2=C(C=CN2)C=C1I",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=CC=C(I)C=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7iodo": {
            "enzyme": "PfTrpB",
            "substrate": "7iodo",
            "substrate-smiles": "C1=CC2=C(C(=C1)I)NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C([O-])C(\C=[NH+]\[C@@H](CC2=CNC3=C(I)C=CC=C23)C([O-])=O)=C(COP([O-])([O-])=O)C=N1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7methyl": {
            "enzyme": "PfTrpB",
            "substrate": "7methyl",
            "substrate-smiles": "CC1=CC=CC2=C1NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "transition-state": ["CC1=C2NC=C(C[C@H](\[NH+]=C\C3=C(COP([O-])([O-])=O)C=NC(C)=C3[O-])C([O-])=O)C2=CC=C1", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        }
    })