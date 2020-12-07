from mhclovac.preprocessing import get_features, transform_ic50_measures
from mhclovac.models import BindingModel
import pandas as pd


TRAINING_SET_SIZE_THRESHOLD = 50
TRAINING_SET_FRACTION = 0.8
MHC_ALLELE = 'HLA-A*02:01'

data = pd.read_csv(f'data/combined_data.tsv', sep='\t', low_memory=False)
