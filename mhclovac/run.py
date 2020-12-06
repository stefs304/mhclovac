from .sequence import validate_sequence, chop_sequence, read_fasta
from .utils import load_model, load_index_data
from .preprocessing import get_features
from .config import Config
import pandas as pd
import sys
import argparse


def predict(sequence: list, mhc: str, sequence_name: str = "unknown") -> pd.DataFrame :
    """
    Predicts binding and epitope score. Returns pandas DataFrame.

    :param sequence: list, list of peptide sequences
    :param mhc: string, MHC allele
    :param sequence_name: sequence_name, defaults to "unknown"
    :return: pandas DataFrame
    """
    model = load_model(mhc)['binding_model']
    data = pd.DataFrame()
    data['sequence'] = sequence
    data['mhc'] = mhc
    data['peptide_length'] = data['sequence'].apply(len)
    data['sequence_name'] = sequence_name
    X = get_features(peptide_list=data['sequence'], index_id_list=Config.INDEX_ID_LIST)
    data['binding_score'] = model.predict(X)
    # sort by combined_score
    data.sort_values(by='binding_score', inplace=True, ascending=False)
    return data


def parse_args(argv):
    description = """
MHCLovac
--------
MHC class I binding and epitope prediction based on modeled physicochemical
properties of peptides. https://gitlab.com/stojanovicbg/mhclovac
Version: 3.3
Author: Stefan Stojanovic
    """
    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-s', '--sequence', type=str, help='Sequence')
    parser.add_argument('-n', '--sequence_name', type=str, help='Sequence name, defaults to Unknown', default='Unknown')
    parser.add_argument('-f', '--fasta', type=str, help='Fasta file')
    parser.add_argument('-m', '--mhc', type=str, help='MHC allele', required=True)
    parser.add_argument('-l', '--peptide_length', type=int, help='Peptide length', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output file name. By default prints to STDOUT')

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
                preds = predict(sequence=peptide_list, mhc=args.mhc, sequence_name=seq_name)
                predictions.append(preds)
            except Exception as e:
                msg = f'Error encountered while processing sequence "{seq_name}": {e}'
                sys.exit(msg)

    if args.sequence:
        try:
            validate_sequence(args.sequence, args.sequence_name, silent=False)
            peptide_list = chop_sequence(args.sequence, args.peptide_length)
            preds = predict(sequence=peptide_list, mhc=args.mhc, sequence_name=args.sequence_name)
            predictions.append(preds)
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
