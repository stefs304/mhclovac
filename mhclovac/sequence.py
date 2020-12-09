import numpy as np
from .utils import pdf


def model_distribution(sequence: str, encoding_scheme: dict, overlap_distance: int = 1,
                       sigma: float = 0.8, n_discrete_points: int = None):
    """
    Models distribution of physicochemical property across a peptide sequence.
    """
    multiplier = 20
    sequence = sequence.upper()
    dist_vector = np.zeros(multiplier*len(sequence)+2*overlap_distance*multiplier)
    for i, A in enumerate(sequence):
        try:
            value = encoding_scheme[A]
        except KeyError:
            msg = 'Unrecognized amino acid: {}'.format(A)
            raise KeyError(msg)
        x = np.linspace(-2.3263, 2.3263, (2*overlap_distance+1)*multiplier)
        A_dist = pdf(x, sigma) * value
        dist_vector[int(i*multiplier):int((i+(2*overlap_distance+1))*multiplier)] += A_dist
    # trim leading and trailing slices
    dist_vector = dist_vector[overlap_distance*multiplier:-overlap_distance*multiplier]
    if n_discrete_points:
        discrete_vector = []
        step = int(len(dist_vector) / n_discrete_points)
        for i in range(0, len(dist_vector), step):
            discrete_vector.append(dist_vector[i])
        return discrete_vector
    return dist_vector


def encode_sequence(sequence: str, encoding_scheme: dict):
    """
    Encodes a peptide sequence with values provided by the encoding table/scheme.
    """
    encoded_sequence = []
    for aa in sequence:
        try:
            value = encoding_scheme.get(aa)
        except Exception as e:
            msg = f'{e}'
            raise KeyError(msg)
        encoded_sequence.append(value)
    return encoded_sequence


def validate_sequence(sequence: str, sequence_name: str = None, silent: bool = True) -> bool:
    """
    Checks if sequence contains only valid amino acid letter codes.

    :param sequence: Amino acid sequence in capital letters
    :param silent: Raises exception if False
    :return: bool, Exception
    """
    valid = ['C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N', 'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M']
    for i, s in enumerate(sequence):
        if s not in valid:
            if not silent:
                msg = f'"{sequence_name}": amino acid "{s}" at position {i} is not valid.'
                raise ValueError(msg)
            return False
    return True


def chop_sequence(sequence: str, peptide_length: int) -> list:
    """
    Chops sequence into N fragments of peptide_length.

    :param sequence: Amino acid sequence
    :param peptide_length: Length of fragments
    :return: list
    """
    if len(sequence) < peptide_length:
        raise ValueError(f'sequence is shorter than peptide_length')
    fragments = []
    for i in range(len(sequence) - peptide_length + 1):
        peptide = sequence[i: i+peptide_length]
        fragments.append(peptide)
    return fragments


def read_fasta(fpath: str) -> (str, str):
    """
    Generator that reads fasta file and yields sequence name and sequence.

    :param fasta_path: string, path to fasta file
    :return: None
    """
    with open(fpath, 'r') as f:
        seq_name = ''
        sequence = ''
        write_flag = 0
        for line in f:
            if line.startswith('>'):
                if write_flag:
                    write_flag = 0
                    yield seq_name, sequence
                seq_name = line[1:].strip()
                sequence = ''
            else:
                sequence += line.strip()
                write_flag = 1
        yield seq_name, sequence


