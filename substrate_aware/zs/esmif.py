"""
A script for generating fasta files for ESMIF.
"""

from __future__ import annotations

import os
from glob import glob
from copy import deepcopy

from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
import pandas as pd
from pathlib import Path
import torch
from tqdm import tqdm

import esm
import esm.inverse_folding

from substrate_aware.preprocess import ZSData
from substrate_aware.util import (
    checkNgen_folder,
    pdb2seq,
    find_missing_str,
    alignmutseq2pdbseq,
    get_chain_structure,
    remove_hetatm,
)


class ESMIFData(ZSData):
    def __init__(
        self,
        input_csv: str,
        in_structure_dir: str = "data/structure",  # data/structure for frompdb or data/structure/docked for docked
        var_col_name: str = "var",
        seq_col_name: str = "seq",
        fit_col_name: str = "fitness",
        seq_dir: str = "data/seq",
        withsub: bool = True,
        esmif_dir: str = "zs/esmif",
        chain_id: str = "A",
    ):

        super().__init__(
            input_csv=input_csv,
            var_col_name=var_col_name,
            seq_col_name=seq_col_name,
            fit_col_name=fit_col_name,
            withsub=withsub,
            seq_dir=seq_dir,
        )

        if self._withsub:
            dets = "docked"
            if in_structure_dir == "data/structure":
                in_structure_dir = "data/structure/docked"
        else:
            dets = "apo"

        self._in_structure_dir = in_structure_dir

        self._esmif_dir = checkNgen_folder(os.path.join(esmif_dir, dets))
        self._esmif_mut_dir = checkNgen_folder(
            os.path.join(self._esmif_dir, "mut_file")
        )

        self._esmif_instruct_dir = checkNgen_folder(
            os.path.join(self._esmif_dir, "instruct_file")
        )
        self._esmif_score_dir = checkNgen_folder(os.path.join(self._esmif_dir, "score"))

        self._chain_id = chain_id

        print(f"Get clean structure {self.esmif_instruct_file}...")
        self._clean_struct()

        print(f"Generating mut fasta file for {self.lib_name}...")
        self._mut_csv2fasta()

        print(f"Scoring mut file for {self.lib_name}...")
        self._score_mut_file()

        score_df = pd.read_csv(self.esmif_output_file)

        sele_cols = [self._var_col_name, self._fit_col_name]

        if "selectivity" in self.df.columns:
            sele_cols.append("selectivity")

        # merge with the original dataframe
        self._esmif_df = pd.merge(
            self.df[sele_cols], score_df, on=self._var_col_name, how="left"
        )

        self._esmif_df.to_csv(self.esmif_output_file, index=False)

    def _clean_struct(self) -> None:
        # TODO direclty using util function or or apo structures from the subfolder
        if self._withsub:
            input_struct = os.path.join(self._in_structure_dir, f"{self.lib_name}.cif")
        else:
            input_struct = os.path.join(
                self._in_structure_dir, f"{self.zs_struct_name}.pdb"
            )

        temp_path = self.esmif_instruct_file.replace(".pdb", "_temp.pdb")

        # get chain structure
        get_chain_structure(
            input_file_path=input_struct,
            output_file_path=temp_path,
            chain_id=self._chain_id,
        )

        remove_hetatm(input_pdb=temp_path, output_pdb=self.esmif_instruct_file)

    def _mut_csv2fasta(self) -> None:
        """
        A function for converting mutation csv to fasta
        """

        # make sure the seq and mut names are in the csv
        for col in [self._var_col_name, self._seq_col_name]:
            if col not in self.df.columns:
                raise ValueError(f"{col} column not found")

        print(f"Writing to {self.esmif_mut_file}...")

        pdb_seq = pdb2seq(self.esmif_instruct_file, self._chain_id)

        # TODO clean up and improve
        if len(self.parent_seq) < len(pdb_seq):
            print("PDB seq is longer than fasta")
            part_before, part_after = find_missing_str(
                longer=pdb_seq, shorter=self.parent_seq
            )
            with open(self.esmif_mut_file, "w") as f:
                for mut, seq in zip(
                    self.df[self._var_col_name].values,
                    self.df[self._seq_col_name].values,
                ):
                    f.write(f">{mut}\n{part_before+seq+part_after}\n")
        elif len(self.parent_seq) == len(pdb_seq):
            print("PDB seq length is equal to fasta")
            with open(self.esmif_mut_file, "w") as f:
                for mut, seq in zip(
                    self.df[self._var_col_name].values,
                    self.df[self._seq_col_name].values,
                ):
                    f.write(f">{mut}\n{seq}\n")
        else:
            print("Fasta seq is longer than PDB")
            index_list, aligned_pdb_seq = alignmutseq2pdbseq(
                mut_seq=self.parent_seq, pdb_seq=pdb_seq
            )

            start_index, end_index = index_list

            # there might be seq with X from pdb
            # Step 1: Find all indices of 'X' in pdb_seq
            x_indices = [i for i, letter in enumerate(aligned_pdb_seq) if letter == "X"]

            with open(self.esmif_mut_file, "w") as f:
                for mut, seq in zip(
                    self.df[self._var_col_name].values,
                    self.df[self._seq_col_name].values,
                ):
                    # Step 2: Modify the original seq by replacing characters at the found indices with 'X'
                    if len(x_indices) > 0 and x_indices[-1] > start_index:
                        start_index = x_indices[-1] + 1
                        # seq_list = list(seq)  # Convert the sequence to a list to allow mutation
                        # for idx in x_indices:
                        #     seq_list[idx] = 'X'

                        # # Convert the modified list back to a string
                        # seq = ''.join(seq_list)
                    f.write(f">{mut}\n{seq[start_index:end_index+1]}\n")

    def _score_mut_file(self) -> None:

        """
        A function for scoring the mutation file
        """
        model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        model = model.eval()

        if torch.cuda.is_available():
            model = model.cuda()
            print("Transferred model to GPU")

        coords, native_seq = esm.inverse_folding.util.load_coords(
            self.esmif_instruct_file, self._chain_id
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
            f"{self.zs_struct_name}.pdb",
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

    @property
    def esmif_df(self) -> pd.DataFrame:
        return self._esmif_df


def run_esmif(
    pattern: str | list = "data/meta/scale2parent/*.csv",
    in_structure_dir: str = "data/structure",
    withsub: bool = False,
    kwargs: dict = {},
):
    """
    Run the esmif scoring for all libraries

    Args:
    - pattern: str | list: the pattern for the input csv files
    """

    if isinstance(pattern, str):
        lib_list = sorted(glob(pattern))
    else:
        lib_list = deepcopy(pattern)

    for lib in lib_list:
        print(f"Running ESMIF for {lib}...")
        ESMIFData(
            input_csv=lib, in_structure_dir=in_structure_dir, withsub=withsub, **kwargs
        )