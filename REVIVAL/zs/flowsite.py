# File to automatically run Flowsite as a ZS
# Author: Lukas Radtke {lradtke@caltech.edu}
# Date: 06/12/24

import numpy as np
import pandas as pd
import os 
import subprocess
import re
import ipdb
import glob
from copy import deepcopy
from tqdm import tqdm
from scipy import stats

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class FlowsiteData(ZSData):
    def __init__(
            self,

    ):
        
        super.__init__()
        


def run_flowsite(pattern: str | list = None, kwargs: dict = {}):

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f'Running Flowsite for {lib}...')
        FlowsiteData(input_csv=lib, **kwargs)