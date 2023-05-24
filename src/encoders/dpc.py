# --------------------------------------------------------------------
# Original code copyright Nicolas Nemeth 2023
# Covered by original MIT license
# --------------------------------------------------------------------


import re
import numpy as np
from typing import List, Tuple


def DPC(fastas: np.ndarray, diPeptides: List[str] = None, seq_range: Tuple[int, int] = None) -> np.ndarray:
    """
        Compute the dipeptide composition.

        Args:
            diPeptides (List[str]): list of dipeptides to compute the composition for
            seq_range (Tuple[int, int]): sequence region to compute dpc for

        Returns:
            n x 400 dimensional np.ndarray where n is the number of protein sequences
    """
    features = list()
    if diPeptides is None:
        diPeptides = [
            aa1 + aa2 for aa1 in "ACDEFGHIKLMNPQRSTVWY" for aa2 in "ACDEFGHIKLMNPQRSTVWY"]
    DPCdict = {dp: 0 for dp in diPeptides}

    if seq_range is not None:
        for idx in range(len(fastas)):
            if len(fastas[idx][1]) > seq_range[1]-1:
                fastas[idx][1] = fastas[idx][1][seq_range[0]:seq_range[1]]

    for _, sequence in fastas:
        sequence = re.sub('-', '', sequence)
        for j in range(len(sequence) - 1):
            if sequence[j]+sequence[j+1] in DPCdict:
                DPCdict[sequence[j]+sequence[j+1]] += 1
        if len(sequence)-1 != 0:
            for dp in DPCdict:
                DPCdict[dp] /= len(sequence)-1
        features.append(list(DPCdict.values()))
    return np.array(features)
