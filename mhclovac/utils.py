

def validate_sequence(peptide):
    # expects uppercase letters
    peptide = peptide.strip()
    if peptide == '':
        return False
    valid = ['C', 'D', 'S', 'Q', 'K', 'I', 'P', 'T', 'F', 'N',
             'G', 'H', 'L', 'R', 'W', 'A', 'V', 'E', 'Y', 'M']
    return not any([p not in valid for p in peptide])



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

