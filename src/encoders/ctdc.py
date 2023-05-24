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


def Count(seq1: str, seq2: str) -> int:
    """
        Compute the total number of all amino acids in seq1 occurring in seq2.

        Args:
            seq1 (str): protein sequence (region) 
            seq2 (str): protein sequence (region)

        Returns:
            int: total count
    """
    sum_ = 0
    for aa in seq1:
        sum_ += seq2.count(aa)
    return sum_


def CTDC(fastas: np.ndarray, seq_range: Tuple[int, int] = None) -> np.ndarray:
    """
        C omposition T ransition D istribution C omposition

        Physicochemical property-based encoding. A method that
        categorizes the 20 primary amino acids into three main classes
        according to seven types of physicochemical properties.

        Layman's terms: percentage of particular amino acid property groups

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
    property = ('hydrophobicity_PRAM900101', 'normwaalsvolume', 'polarity',
                'polarizability', 'charge', 'secondarystruct', 'solventaccess')

    if seq_range is not None:
        for idx in range(len(fastas)):
            if len(fastas[idx][1]) > seq_range[1]-1:
                fastas[idx][1] = fastas[idx][1][seq_range[0]:seq_range[1]]

    features = list()
    for _, sequence in fastas:
        code = list()
        sequence = re.sub('-', '', sequence)
        for p in property:
            c1 = Count(group1[p], sequence) / len(sequence)
            c2 = Count(group2[p], sequence) / len(sequence)
            c3 = Count(group3[p], sequence) / len(sequence)
            code.append(c1)
            code.append(c2)
            code.append(c3)
        features.append(code)
    return np.array(features)
