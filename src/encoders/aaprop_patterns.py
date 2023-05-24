# --------------------------------------------------------------------
# Original code copyright Nicolas Nemeth 2023
# Covered by original MIT license
# --------------------------------------------------------------------


import numpy as np
from typing import List, Tuple

# Pattern Data
POLAR = "NQST"
AROMATIC = "VWY"
ACIDIC = "DE"
ALKALINE = "KLR"
IONISABLE = "CY"
HYDROPHOBIC_2 = "VMLAIG"
HYDROPHILIC_2 = "SFTNKYEQCWPHDRU"

HYDROPHOBIC_HYDROPHILIC_2NDALPHABET = [[HYDROPHILIC_2, HYDROPHILIC_2],
                                       [HYDROPHILIC_2, HYDROPHILIC_2, HYDROPHILIC_2],
                                       [HYDROPHILIC_2, HYDROPHILIC_2, HYDROPHOBIC_2],
                                       [HYDROPHILIC_2, HYDROPHOBIC_2, HYDROPHILIC_2],
                                       [HYDROPHILIC_2, HYDROPHOBIC_2, HYDROPHOBIC_2],
                                       [HYDROPHOBIC_2, HYDROPHILIC_2, HYDROPHILIC_2],
                                       [HYDROPHOBIC_2, HYDROPHILIC_2, HYDROPHOBIC_2],
                                       [HYDROPHOBIC_2, HYDROPHOBIC_2, HYDROPHILIC_2]]

AA_PROPERTY_2NDALPHABET = [[ACIDIC], [ACIDIC, HYDROPHOBIC_2], [ALKALINE], [ALKALINE, ALKALINE],
                           [ALKALINE, HYDROPHILIC_2], [
                               ALKALINE, HYDROPHOBIC_2],
                           [ALKALINE, HYDROPHOBIC_2, HYDROPHOBIC_2], [
                               ALKALINE, HYDROPHOBIC_2, POLAR],
                           [ALKALINE, POLAR], [ALKALINE, POLAR, POLAR], [
                               AROMATIC], [HYDROPHILIC_2],
                           [HYDROPHILIC_2, ALKALINE], [HYDROPHILIC_2,
                                                       HYDROPHOBIC_2], [HYDROPHILIC_2, POLAR],
                           [HYDROPHOBIC_2], [HYDROPHOBIC_2, ACIDIC], [
                               HYDROPHOBIC_2, ALKALINE],
                           [HYDROPHOBIC_2, ALKALINE, HYDROPHOBIC_2], [
                               HYDROPHOBIC_2, ALKALINE, POLAR],
                           [HYDROPHOBIC_2, HYDROPHILIC_2], [
                               HYDROPHOBIC_2, HYDROPHOBIC_2],
                           [HYDROPHOBIC_2, HYDROPHOBIC_2, HYDROPHOBIC_2], [
                               HYDROPHOBIC_2, HYDROPHOBIC_2, POLAR],
                           [HYDROPHOBIC_2, IONISABLE], [HYDROPHOBIC_2, POLAR], [
                               HYDROPHOBIC_2, POLAR, HYDROPHOBIC_2],
                           [HYDROPHOBIC_2, POLAR, POLAR], [IONISABLE], [
                               IONISABLE, POLAR], [POLAR], [POLAR, ACIDIC],
                           [POLAR, ALKALINE], [POLAR, ALKALINE, HYDROPHOBIC_2], [
                               POLAR, ALKALINE, POLAR],
                           [POLAR, HYDROPHILIC_2], [
                               POLAR, HYDROPHILIC_2, HYDROPHOBIC_2],
                           [POLAR, HYDROPHOBIC_2, HYDROPHOBIC_2], [
                               POLAR, HYDROPHOBIC_2, POLAR],
                           [POLAR, POLAR], [POLAR, POLAR, ALKALINE], [
                               POLAR, POLAR, HYDROPHOBIC_2],
                           [POLAR, POLAR, POLAR]]

PATTERN_LABELS = ["HYDROPHILIC_2, HYDROPHILIC_2",
                  "HYDROPHILIC_2, HYDROPHILIC_2, HYDROPHILIC_2",
                  "HYDROPHILIC_2, HYDROPHILIC_2, HYDROPHOBIC_2",
                  "HYDROPHILIC_2, HYDROPHOBIC_2, HYDROPHILIC_2",
                  "HYDROPHILIC_2, HYDROPHOBIC_2, HYDROPHOBIC_2",
                  "HYDROPHOBIC_2, HYDROPHILIC_2, HYDROPHILIC_2",
                  "HYDROPHOBIC_2, HYDROPHILIC_2, HYDROPHOBIC_2",
                  "HYDROPHOBIC_2, HYDROPHOBIC_2, HYDROPHILIC_2"] + ["ACIDIC", "ACIDIC, HYDROPHOBIC_2",
                                                                    "ALKALINE", "ALKALINE, ALKALINE",
                                                                    "ALKALINE, HYDROPHILIC_2", "ALKALINE, HYDROPHOBIC_2",
                                                                    "ALKALINE, HYDROPHOBIC_2, HYDROPHOBIC_2", "ALKALINE, HYDROPHOBIC_2, POLAR",
                                                                    "ALKALINE, POLAR", "ALKALINE, POLAR, POLAR", "AROMATIC", "HYDROPHILIC_2",
                                                                    "HYDROPHILIC_2, ALKALINE", "HYDROPHILIC_2, HYDROPHOBIC_2", "HYDROPHILIC_2, POLAR",
                                                                    "HYDROPHOBIC_2", "HYDROPHOBIC_2, ACIDIC", "HYDROPHOBIC_2, ALKALINE",
                                                                    "HYDROPHOBIC_2, ALKALINE, HYDROPHOBIC_2", "HYDROPHOBIC_2, ALKALINE, POLAR",
                                                                    "HYDROPHOBIC_2, HYDROPHILIC_2", "HYDROPHOBIC_2, HYDROPHOBIC_2",
                                                                    "HYDROPHOBIC_2, HYDROPHOBIC_2, HYDROPHOBIC_2", "HYDROPHOBIC_2, HYDROPHOBIC_2, POLAR",
                                                                    "HYDROPHOBIC_2, IONISABLE", "HYDROPHOBIC_2, POLAR", "HYDROPHOBIC_2, POLAR, HYDROPHOBIC_2",
                                                                    "HYDROPHOBIC_2, POLAR, POLAR", "IONISABLE", "IONISABLE, POLAR", "POLAR", "POLAR, ACIDIC",
                                                                    "POLAR, ALKALINE", "POLAR, ALKALINE, HYDROPHOBIC_2", "POLAR, ALKALINE, POLAR",
                                                                    "POLAR, HYDROPHILIC_2", "POLAR, HYDROPHILIC_2, HYDROPHOBIC_2",
                                                                    "POLAR, HYDROPHOBIC_2, HYDROPHOBIC_2", "POLAR, HYDROPHOBIC_2, POLAR",
                                                                    "POLAR, POLAR", "POLAR, POLAR, ALKALINE", "POLAR, POLAR, HYDROPHOBIC_2",
                                                                    "POLAR, POLAR, POLAR"]


SELECTED_PATTERNS = HYDROPHOBIC_HYDROPHILIC_2NDALPHABET + AA_PROPERTY_2NDALPHABET


def AaPropPatterns(fastas: np.ndarray, patterns: List[List[str]] = SELECTED_PATTERNS,
                   seq_range: Tuple[int, int] = None) -> np.ndarray:
    """
        Input numpy array of [protein name, sequence] entries.
        Outputs matrix where each row is a vector containing the frequency of
        each pattern in the protein sequence.

        Args:
            fastas: numpy array of [protein name, sequence] entries.
            patterns: list of patterns to search for in the sequence.
            seq_range: tuple of (start, end) indices to extract from the sequence.

        Returns:
            results: numpy array of [protein name, vector of pattern frequencies] entries.
    """
    if seq_range is not None:
        sequences = [seq[seq_range[0]:seq_range[1]] if len(
            seq) > seq_range[1] else seq for seq in fastas[:, 1]]
    else:
        sequences = fastas[:, 1]

    results = [[compute_pattern_freq(seq, pattern)
                for pattern in patterns] for seq in sequences]
    return np.array(results)


def compute_pattern_freq(seq: str, pattern: str) -> float:
    """
        Compute the frequency of a given pattern in a protein sequence.

        The frequency is computed as the count of non-overlapping instances of the pattern
        divided by the maximum possible count given the length of the sequence and pattern.

        Args:
            seq (str): A protein sequence, represented as a string of amino acids.
            pattern (str): A pattern to be searched for in the sequence. This is a list of amino acid groups,
                     where each group is a string of one or more possible amino acids at that position in the pattern.

        Returns:
            float: The frequency of the pattern in the sequence, computed as described above.
            If the pattern is longer than the sequence, returns 0.
    """
    max_p_count = len(seq) - len(pattern) + 1
    if max_p_count > 0:
        p_count = sum(all(aa in group for aa, group in zip(
            seq[pos:pos+len(pattern)], pattern)) for pos in range(max_p_count))
        return p_count / max_p_count
    return 0
