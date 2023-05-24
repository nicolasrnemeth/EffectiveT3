# --------------------------------------------------------------------
# Original code copyright Zhen Chen, Xuhan Liu, Pei Zhao, Chen Li, Yanan Wang, Fuyi Li, Tatsuya Akutsu, Chris Bain, Robin B Gasser, Junzhou Li, Zuoren Yang, Xin Gao, Lukasz Kurgan, Jiangning Song 2022
# Covered by original MIT license
##
# Modifications copyright Nicolas Nemeth 2023
# Modifications licensed under the MIT License
# --------------------------------------------------------------------


import re
import numpy as np
from typing import Tuple


def CTDT(fastas: np.ndarray, seq_range: Tuple[int, int] = None) -> np.ndarray:
    """
        C omposition T ransition D istribution T ransition

        Physicochemical property-based encoding. A method that
        categorizes the 20 primary amino acids into three main classes
        according to seven types of physicochemical properties.

        Layman's terms: percentage of mutual conversion in amino acid properties

        Args:
            fastas (np.ndarray): array containing 2-sized list -> protein identifier | protein sequence
            seq_range (Tuple[int, int]): sequence region to compute the encoding for

        Returns:
            n x 21 dimensional np.ndarray where n is the number of protein sequences
    """

    group1 = {
        'hydrophobicity_PRAM900101': 'RKEDQN',
        'normwaalsvolume': 'GASTPDC',
        'polarity':        'LIFWCMVY',
        'polarizability':  'GASDT',
        'charge':          'KR',
        'secondarystruct': 'EALMQKRH',
        'solventaccess':   'ALFCGIVW'
    }
    group2 = {
        'hydrophobicity_PRAM900101': 'GASTPHY',
        'normwaalsvolume': 'NVEQIL',
        'polarity':        'PATGS',
        'polarizability':  'CPNVEQIL',
        'charge':          'ANCQGHILMFPSTWYV',
        'secondarystruct': 'VIYCWFT',
        'solventaccess':   'RKQEND'
    }
    group3 = {
        'hydrophobicity_PRAM900101': 'CLVIMFW',
        'normwaalsvolume': 'MHKFRYW',
        'polarity':        'HQRKNED',
        'polarizability':  'KMHFRYW',
        'charge':          'DE',
        'secondarystruct': 'GNPSD',
        'solventaccess':   'MSPTHY'
    }

    property = (
        'hydrophobicity_PRAM900101', 'normwaalsvolume',
        'polarity', 'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    if seq_range is not None:
        for idx in range(len(fastas)):
            if len(fastas[idx][1]) > seq_range[1]-1:
                fastas[idx][1] = fastas[idx][1][seq_range[0]:seq_range[1]]

    features = list()
    for _, sequence in fastas:
        code = list()
        sequence = re.sub('-', '', sequence)
        aaPair = [sequence[j:j + 2] for j in range(len(sequence) - 1)]
        for p in property:
            c1221, c1331, c2332 = 0, 0, 0
            for pair in aaPair:
                if (pair[0] in group1[p] and pair[1] in group2[p]) or (pair[0] in group2[p] and pair[1] in group1[p]):
                    c1221 += 1
                    continue
                if (pair[0] in group1[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group1[p]):
                    c1331 += 1
                    continue
                if (pair[0] in group2[p] and pair[1] in group3[p]) or (pair[0] in group3[p] and pair[1] in group2[p]):
                    c2332 += 1
            code.append(c1221/len(aaPair))
            code.append(c1331/len(aaPair))
            code.append(c2332/len(aaPair))
        features.append(code)
    return np.array(features)
