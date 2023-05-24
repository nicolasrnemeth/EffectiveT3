# --------------------------------------------------------------------
# Original code copyright Nicolas Nemeth 2023
# Covered by original MIT license
# --------------------------------------------------------------------


import re
import numpy as np
from typing import Tuple


def AAC(fastas: np.ndarray, AAs: str = None, seq_range: Tuple[int, int] = None) -> np.ndarray:
    """
        Compute the amino acid composition for all 
        amino acids provided with the parameter ´AAs´.

        Args:
            fastas (np.ndarray): array of 2-sized lists -> protein identifier | protein sequence
            AAs (str): amino acids to compute the amino acid composition for
            seq_range (Tuple[int, int]): sequence region to compute the amino acid composition for

        Returns:
            array with encoded amino acid composition of the protein sequence
    """
    if seq_range is not None:
        for idx in range(len(fastas)):
            if len(fastas[idx][1]) > seq_range[1]-1:
                fastas[idx][1] = fastas[idx][1][seq_range[0]:seq_range[1]]
    if AAs is None:
        AAs = "ACDEFGHIKLMNPQRSTVWY"

    features = list()
    for _, sequence in fastas:
        sequence = re.sub('-', '', sequence)
        features.append([sequence.count(aa) / len(sequence) for aa in AAs])
    return np.array(features)
