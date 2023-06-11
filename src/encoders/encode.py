import numpy as np
from typing import Tuple
from .aac import AAC
from .ctdc import CTDC
from .ctdt import CTDT
from .dpc import DPC
from .aaprop_patterns import AaPropPatterns

# Dipeptides used for DPC computation as obtained by feature selection using
# shap-values of features obtained from the test set of trained light gradient boosting model
# Split in 3 groups to maintain the same order as obtained by feature selection
# otherwise the same model with that performs so well on multiple test sets
# cannot be obtained again with the exact same optimized weight (parameters)
DPC_FEATURE_SELECTION_1 = [
    'SS' 'KR' 'PS' 'WV' 'LH' 'TP' 'FY' 'WG' 'WA' 'NS' 'SN' 'SP' 'PI' 'WL'
    'PP' 'NT' 'AR' 'PT' 'NN' 'FA' 'EW' 'IW' 'VE' 'VV' 'VI' 'QS' 'VL' 'QP'
    'IL' 'TS' 'QN' 'ER' 'LW' 'WR' 'II' 'SQ' 'ST' 'RD' 'SC' 'GF' 'TQ' 'LM'
    'HS' 'WD' 'SG'
]
DPC_FEATURE_SELECTION_2 = [
    'HT' 'DW' 'TT' 'QM' 'AP' 'QH' 'TN' 'LV' 'FE' 'LA'
    'AW' 'PW' 'SH' 'VD' 'RG' 'WQ' 'QT' 'DE' 'KW' 'DF' 'NH'
]
DPC_FEATURE_SELECTION_3 = [
    'EV' 'FG' 'VM' 'RS' 'LL' 'AL' 'DY' 'AI' 'IV' 'FI'
    'IA' 'YS' 'PA' 'DI' 'IN' 'TL' 'NP'
]
# AaProp feature selection using feature importances of trained light gradient boosting model
POLAR = "NQST"
# Second set of dipeptides for DPC as obtained by feature selection using feature
# importances of trained light gradient boosting model
# -> hardcoded: see function ´src/encoders/ctdc.py´
# uncomment code to obtain full ctdc encoding


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
             -> returns n x 85 dimensional np.ndarray, 85 are the shortlisted
             features as obtained by the feature selection process using shap-values,
             on different test sets
    """

    # Only perform sequence region extraction once to improve computation time
    if seq_range is not None:
        for idx in range(len(fastas)):
            if len(fastas[idx][1]) > seq_range[1]-1:
                fastas[idx][1] = fastas[idx][1][seq_range[0]:seq_range[1]]

    # Sequence-based features
    dpc1 = DPC(fastas, seq_range=None, diPeptides=DPC_FEATURE_SELECTION_1)
    dpc2 = DPC(fastas, seq_range=None, diPeptides=DPC_FEATURE_SELECTION_2)
    dpc3 = DPC(fastas, seq_range=None, diPeptides=DPC_FEATURE_SELECTION_3)
    # Amino acid property (patterns) based features
    aaprop = AaPropPatterns(fastas, seq_range=None, patterns=[POLAR])
    # Feature selection is hardcoded inside the function ´src/encoders/ctdc.py´
    ctdc = CTDC(fastas, seq_range=None)
    # Combine all features
    features = np.hstack((dpc1, aaprop, dpc2, ctdc, dpc3))
    # names, encodings
    return fastas[:, 0], features


# Obtain the full encoding including all features contained in the encoders module
# copy the function body of ´full_encode´ and paste it into the ´encode´ function
# if you want to train a model using all features to later use it for prediction
# Note the model with less features was shown to perform better

# NOTE REMOVE THE COMMENTS IN THE FUNCTION ´src/encoders/ctdc.py´ TO OBTAIN THE FULL CTDC ENCODING
# FOR REASON OF EFFICIENCY THE FEATURE SELECTED ENCODING OF CTDC IS HARDCODED IN THE FUNCTION

# def full_encode(fastas: np.ndarray, seq_range: Tuple[int, int] = None) -> Tuple[np.ndarray, np.ndarray]:
#     """
#         Encode protein sequences using sequence- and amino acid property-based features.

#         Args:
#             fastas (np.ndarray): array containing the protein sequences
#             seq_range (Tuple[int, int]): sequence range to use for prediction (defaults to full-length)
#             model_ (str): single model combination to use for prediction
#             process_id (str): id of the process, in case this task is executed in parallel
#             result_dict (dict): dictionary to add the results to

#         Returns:
#             Tuple[np.ndarray, np.ndarray]: containing the protein identifiers
#             and the encoded features of all input protein sequences
#     """
#     # Sequence-based features
#     aac = AAC(fastas, seq_range=seq_range)
#     dpc = DPC(fastas, seq_range=seq_range)
#     # Amino acid property (patterns) based features
#     aaprop = AaPropPatterns(fastas, seq_range=seq_range)
#     ctdc = CTDC(fastas, seq_range=seq_range)
#     ctdt = CTDT(fastas, seq_range=seq_range)
#     # Combine all features
#     features = np.hstack((aac, aaprop, ctdc, ctdt, dpc))
#     # names, encodings
#     return fastas[:, 0], features
