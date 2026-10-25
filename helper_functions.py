import numpy as np

def jukes_cantor(reference_sequence: str, distant_sequence: str) -> float:
    """The Jukes-Cantor correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Jukes-Cantor corrected genetic distance using Hamming distance.
        For example 1.163.

    """
    mismatches = 0
    valid_positions = 0

    for ref, seq in zip(reference_sequence, distant_sequence):
        if ref not in "ACGT" or seq not in "ACGT":
            continue

        valid_positions += 1

        if ref != seq:
            mismatches += 1

    if valid_positions == 0:
        return 0.0

    p = mismatches / valid_positions
    if p >= 0.75:
            return float("inf")
    
    jc_distance = -3/4 * np.log(1 - 4/3 * p)
    answer = jc_distance * valid_positions

    return answer


def kimura_two_parameter(reference_sequence: str, distant_sequence: str) -> float:
    """The Kimura Two Parameter correction for estimating genetic distances
    calculated with Hamming distance.
    Should return genetic distance with the same unit as if not corrected.

    Parameters
    ----------
    reference_sequence: str
        A string of nucleotides in a sequence used as a reference
        in an alignment with other (e.g. AGGT-GA)
    distant_sequence: str
        A string of nucleotides in a sequence after the alignment
        with a reference (e.g. AGC-AGA)

    Returns
    -------
    float
        The Kimura corrected genetic distance using Hamming distance.
        For example 1.196.

    """
    transitions = 0
    transversions = 0
    valid_positions = 0

    for ref, seq in zip(reference_sequence, distant_sequence):
        if ref not in "ACTG" or seq not in "ACTG":
            continue

        valid_positions += 1

        if ref != seq:
            if (ref in "AG" and seq in "AG") or (ref in "CT" and seq in "CT"):
                transitions += 1
            else:
                transversions += 1

    if valid_positions == 0:
        return 0.0

    P = transitions / valid_positions
    Q = transversions / valid_positions

    if 1 - 2 * P - Q <= 0 or 1 - 2 * Q <= 0:
        return float("inf")

    k2p = -0.5 * np.log(1 - 2 * P - Q) - 0.25 * np.log(1 - 2 * Q)

    answer = k2p * valid_positions
    return answer
            
