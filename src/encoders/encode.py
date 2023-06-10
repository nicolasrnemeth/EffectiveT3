import numpy as np
from typing import Tuple
from .aac import AAC
from .ctdc import CTDC
from .ctdt import CTDT
from .dpc import DPC
from .aaprop_patterns import AaPropPatterns

# # DPC as obtained by feature selection using feature importances of trained light gradient boosting model
# DPC_FEATURE_SELECTION = [
#     "LW", "SP", "YF", "VI", "SS", "FA", "NN", "SH", "WL", "NT",
#     "RD", "SN", "WV", "SQ", "EV", "KR", "ST", "TP", "GY", "PS", "PT"
# ]
# # AaProp feature selection using feature importances of trained light gradient boosting model
# POLAR = "NQST"


def encode(fastas: np.ndarray, seq_range: Tuple[int, int] = None) -> Tuple[np.ndarray, np.ndarray]:
    """
        Encode protein sequences using sequence- and amino acid property-based features.

        Args:
            fastas (np.ndarray): array containing the protein sequences
            seq_range (Tuple[int, int]): sequence range to use for prediction (defaults to full-length)
            model_ (str): single model combination to use for prediction
            process_id (str): id of the process, in case this task is executed in parallel
            result_dict (dict): dictionary to add the results to

        Returns:
            Tuple[np.ndarray, np.ndarray]: containing the protein identifiers
            and the encoded features of all input protein sequences
    """
    # Sequence-based features
    aac = AAC(fastas, seq_range=seq_range)
    dpc = DPC(fastas, seq_range=seq_range)
    # Amino acid property (patterns) based features
    aaprop = AaPropPatterns(fastas, seq_range=seq_range)
    ctdc = CTDC(fastas, seq_range=seq_range)
    ctdt = CTDT(fastas, seq_range=seq_range)
    # Combine all features
    features = np.hstack((aac, aaprop, ctdc, ctdt, dpc))
    # names, encodings
    return fastas[:, 0], features
