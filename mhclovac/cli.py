
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', help='Input fasta file.')
    parser.add_argument('--peptides', help='Peptides; comma separated.')
    parser.add_argument('--mhc', help='List of MHC alleles; comma separated.')

