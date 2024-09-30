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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-6bromo": {
            "enzyme": "PfTrpB",
            "substrate": "6bromo",
            "substrate-smiles": "C1=CC(=CC2=C1C=CN2)Br",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        "PfTrpB-7bromo": {
            "enzyme": "PfTrpB",
            "substrate": "7bromo",
            "substrate-smiles": "C1=CC2=C(C(=C1)Br)NC=C2",
            "cofactor": ["PLP", "Na+"],
            "cofactor-smiles": ["O=Cc1c(O)c(C)ncc1COP(O)(O)=O", "[Na+]"],
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
            "positions": {1: 165, 2: 183, 3: 301},
            "AAs": {1: "I", 2: "I", 3: "Y"}, # I165, I183, and Y301
            "family": "TrpB",
            "project": "multi-substrate",
        },
        # below are SSMuLA datasets
        "DHFR": {
            "positions": {1: 26, 2: 27, 3: 28},
            "codons": {1: "GCC", 2: "GAT", 3: "CTC"},
            "AAs": {1: "A", 2: "D", 3: "L"},
            "type": "Enzyme activity",
        },
        "ParD2": {
            "positions": {1: 61, 2: 64, 3: 80},
            "codons": {1: "", 2: "", 3: ""},
            "AAs": {1: "I", 2: "L", 3: "K"},
            "type": "Binding",
        },
        "ParD3": {
            "positions": {1: 61, 2: 64, 3: 80},
            "codons": {1: "", 2: "", 3: ""},
            "AAs": {1: "D", 2: "K", 3: "E"},
            "type": "Binding",
        },
        "GB1": {
            "positions": {1: 39, 2: 40, 3: 41, 4: 54},
            "codons": {1: "", 2: "", 3: "", 4: ""},
            "AAs": {1: "V", 2: "D", 3: "G", 4: "V"},
            "type": "Binding",
        },
        "T7": {
            "positions": {1: 748, 2: 756, 3: 758},
            "codons": {1: "", 2: "", 3: ""},
            "AAs": {1: "N", 2: "R", 3: "Q"},
            "type": "Enzyme activity",
            "pdb_resrange": [71, 883],
        },
        "TEV": {
            "positions": {1: 146, 2: 148, 3: 167, 4: 170},
            "codons": {1: "", 2: "", 3: "", 4: ""},
            "AAs": {1: "T", 2: "D", 3: "H", 4: "S"},
            "type": "Enzyme activity",
            "pdb_resrange": [1, 218],
        },
        "TrpB3A": {
            "positions": {1: 104, 2: 105, 3: 106},
            "codons": {1: "GCT", 2: "GAA", 3: "ACG"},
            "AAs": {1: "A", 2: "E", 3: "T"},
            "type": "Enzyme activity",
        },
        "TrpB3B": {
            "positions": {1: 105, 2: 106, 3: 107},
            "codons": {1: "GAA", 2: "ACG", 3: "GGT"},
            "AAs": {1: "E", 2: "T", 3: "G"},
            "type": "Enzyme activity",
        },
        "TrpB3C": {
            "positions": {1: 106, 2: 107, 3: 108},
            "codons": {1: "ACG", 2: "GGT", 3: "GCT"},
            "AAs": {1: "T", 2: "G", 3: "A"},
            "type": "Enzyme activity",
        },
        "TrpB3D": {
            "positions": {1: 117, 2: 118, 3: 119},
            "codons": {1: "ACC", 2: "GCA", 3: "GCA"},
            "AAs": {1: "T", 2: "A", 3: "A"},
            "type": "Enzyme activity",
        },
        "TrpB3E": {
            "positions": {1: 184, 2: 185, 3: 186},
            "codons": {1: "TTC", 2: "GGC", 3: "TCT"},
            "AAs": {1: "F", 2: "G", 3: "S"},
            "type": "Enzyme activity",
        },
        "TrpB3F": {
            "positions": {1: 162, 2: 166, 3: 301},
            "codons": {1: "CTG", 2: "ATT", 3: "TAC"},
            "AAs": {1: "L", 2: "I", 3: "Y"},
            "type": "Enzyme activity",
        },
        "TrpB3G": {
            "positions": {1: 227, 2: 228, 3: 301},
            "codons": {1: "GTG", 2: "AGC", 3: "TAC"},
            "AAs": {1: "V", 2: "S", 3: "Y"},
            "type": "Enzyme activity",
        },
        "TrpB3H": {
            "positions": {1: 228, 2: 230, 3: 231},
            "codons": {1: "AGC", 2: "GGT", 3: "TCT"},
            "AAs": {1: "S", 2: "G", 3: "S"},
            "type": "Enzyme activity",
        },
        "TrpB3I": {
            "positions": {1: 182, 2: 183, 3: 184},
            "codons": {1: "TAC", 2: "GTG", 3: "TTC"},
            "AAs": {1: "Y", 2: "V", 3: "F"},
            "type": "Enzyme activity",
        },
        "TrpB4": {
            "positions": {1: 183, 2: 184, 3: 227, 4: 228},
            "codons": {1: "GTG", 2: "TTC", 3: "GTG", 4: "AGC"},
            "AAs": {1: "V", 2: "F", 3: "V", 4: "S"},
            "type": "Enzyme activity",
        },
    })