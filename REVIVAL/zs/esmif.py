"""
A script for generating fasta files for ESMIF.
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy

from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
from tqdm import tqdm

import esm
import esm.inverse_folding

from REVIVAL.preprocess import ZSData
from REVIVAL.util import checkNgen_folder


class ESMIFData(ZSData):
    def __init__(
        self,
        input_csv: str,
        scale_fit: str = "parent",
        combo_col_name: str = "AAs",
        var_col_name: str = "var",
        mut_col_name: str = "mut",
        pos_col_name: str = "pos",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        structure_dir: str = "data/structure",
        zs_dir: str = "zs",
        esmif_dir: str = "esmif",
        esmif_mut_dir: str = "mut_file",
        esmif_instruct_dir: str = "input_structure",
        esmif_score_dir: str = "output",
        chain_id: str = "A",
    ):

        super().__init__(
            input_csv=input_csv,
            scale_fit=scale_fit,
            combo_col_name=combo_col_name,
            var_col_name=var_col_name,
            mut_col_name=mut_col_name,
            pos_col_name=pos_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            seq_dir=seq_dir,
            structure_dir=structure_dir,
            zs_dir=zs_dir,
        )

        self._esmif_dir = checkNgen_folder(os.path.join(self._zs_dir, esmif_dir))
        self._esmif_mut_dir = checkNgen_folder(
            os.path.join(self._esmif_dir, esmif_mut_dir)
        )
        self._esmif_instruct_dir = checkNgen_folder(
            os.path.join(self._esmif_dir, esmif_instruct_dir)
        )
        self._esmif_score_dir = checkNgen_folder(
            os.path.join(self._esmif_dir, esmif_score_dir)
        )
        self._chain_id = chain_id

        print(f"Generating mut fasta file for {self.lib_name}...")
        self._mut_csv2fasta()
        print(f"Scoring mut file for {self.lib_name}...")
        self._score_mut_file()

    def _mut_csv2fasta(self) -> None:
        """
        A function for converting mutation csv to fasta
        """

        # make sure the seq and mut names are in the csv
        for col in [self._var_col_name, self._seq_col_name]:
            if col not in self.df.columns:
                raise ValueError(f"{col} column not found")

        print(f"Writing to {self.esmif_mut_file}...")

        # TODO: dealt with seq not same length
        with open(self.esmif_mut_file, "w") as f:
            for mut, seq in zip(
                self.df[self._var_col_name].values, self.df[self._seq_col_name].values
            ):
                f.write(f">{mut}\n{seq}\n")

    def _score_mut_file(self) -> None:

        """
        A function for scoring the mutation file
        """
        model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        model = model.eval()

        if torch.cuda.is_available():
            model = model.cuda()
            print("Transferred model to GPU")
        print(f"in esmifdata self._structure_dir: {self._structure_dir}")
        coords, native_seq = esm.inverse_folding.util.load_coords(
            self.structure_file, self._chain_id
        )
        print(f"Native sequence loaded from structure file:\n{native_seq}\n")

        ll, _ = esm.inverse_folding.util.score_sequence(
            model, alphabet, coords, native_seq
        )
        print("Native sequence")
        print(f"Log likelihood: {ll:.2f}")
        print(f"Perplexity: {np.exp(-ll):.2f}")

        print("\nScoring variant sequences from sequence file..\n")

        infile = FastaFile()
        infile.read(self.esmif_mut_file)
        seqs = get_sequences(infile)
        Path(self.esmif_output_file).parent.mkdir(parents=True, exist_ok=True)
        with open(self.esmif_output_file, "w") as fout:
            fout.write(f"{self._var_col_name},esmif_score\n")
            for header, seq in tqdm(seqs.items()):
                ll, _ = esm.inverse_folding.util.score_sequence(
                    model, alphabet, coords, str(seq)
                )
                fout.write(header + "," + str(ll) + "\n")
        print(f"Results saved to {self.esmif_output_file}")

    @property
    def esmif_instruct_file(self) -> str:
        """
        PDB file path to the esmif pdb file
        """
        return os.path.join(
            self._esmif_instruct_dir,
            f"{self.lib_name}.{os.path.splitext(self.structure_file)[-1]}",
        )

    @property
    def esmif_mut_file(self) -> str:
        """
        Mutation file path to the esmif mutation file
        """
        return os.path.join(self._esmif_mut_dir, f"{self.lib_name}.fasta")

    @property
    def esmif_output_file(self) -> str:
        """
        Score file path to the esmif score file
        """
        return os.path.join(self._esmif_score_dir, f"{self.lib_name}.csv")


def run_esmif(pattern: str | list = "data/meta/scale2parent/*.csv", kwargs: dict = {}):
    """
    Run the esmif scoring for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    """

    if isinstance(pattern, str):
        lib_list = glob(pattern)
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running ESMIF for {lib}...")
        ESMIFData(input_csv=lib, **kwargs)