from .sequence import validate_sequence, chop_sequence, read_fasta
from .utils import load_model
from .preprocessing import get_features
from .config import Config
import pandas as pd
import sys
import argparse


def predict(peptides: list, mhc_allele: str, sequence_name: str = "unknown", sort: bool = False, n_cpu=1) -> pd.DataFrame:
    """
    Predicts binding score (affinity). Returns pandas DataFrame.
    """
    model = load_model(mhc_allele)
    data = pd.DataFrame()
    data['peptide'] = peptides
    data['mhc_allele'] = mhc_allele
    data['peptide_length'] = data['peptide'].apply(len)
    data['sequence_name'] = sequence_name
    x = get_features(peptide_list=data['peptide'], index_id_list=Config.INDEX_ID_LIST, n_cpu=n_cpu)
    data['binding_score'] = model.predict(x)
    # sort values optionally
    if sort:
        data.sort_values(by='binding_score', inplace=True, ascending=False)
    return data


def parse_args(argv):
    description = """
MHCLovac
--------
MHC binding prediction based on modeled physicochemical properties of peptides. 
https://github.com/stefs304/mhclovac
Version: 4.0
Author: Stefan Stojanovic
"""
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--sequence', type=str, help='Sequence')
    parser.add_argument('-n', '--sequence_name', type=str, help='Sequence name', default='Unknown')
    parser.add_argument('-f', '--fasta', type=str, help='Fasta file')
    parser.add_argument('-m', '--mhc', type=str, help='MHC allele', required=True)
    parser.add_argument('-l', '--peptide_length', type=int, help='Peptide length', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name. By default, output is printed to STDOUT')
    parser.add_argument('--sort', action='store_true', help='Sort output based on prediction score')
    parser.add_argument('--n_cpu', type=int, help='Number of CPU cores to use', default=1)

    return parser.parse_args(argv)


def run():
    args = parse_args(sys.argv[1:])

    if not any([args.sequence, args.fasta]):
        msg = f'Sequence or fasta file are required.\n'
        sys.exit(msg)

    predictions = []
    if args.fasta:
        for seq_name, sequence in read_fasta(args.fasta):
            try:
                validate_sequence(sequence, seq_name, silent=False)
                peptide_list = chop_sequence(sequence, args.peptide_length)
                p = predict(peptides=peptide_list, mhc_allele=args.mhc, sequence_name=seq_name, sort=args.sort, n_cpu=args.n_cpu)
                predictions.append(p)
            except Exception as e:
                msg = f'Error encountered while processing sequence "{seq_name}": {e}'
                sys.exit(msg)

    if args.sequence:
        try:
            validate_sequence(args.sequence, args.sequence_name, silent=False)
            peptide_list = chop_sequence(args.sequence, args.peptide_length)
            p = predict(peptides=peptide_list, mhc_allele=args.mhc, sequence_name=args.sequence_name, sort=args.sort, n_cpu=args.n_cpu)
            predictions.append(p)
        except Exception as e:
            msg = f'Error encountered while processing sequence "{args.sequence_name}": {e}'
            sys.exit(msg)

    try:
        output = pd.concat(predictions)
    except Exception as e:
        sys.exit(e)

    if args.output:
        output.to_csv(args.output, sep='\t', index=False)
    else:
        sys.stdout.write(output.to_string(index=False) + '\n')

    sys.exit()


if __name__ == '__main__':
    run()
