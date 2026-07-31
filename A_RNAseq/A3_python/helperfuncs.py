import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from scipy import stats

from statsmodels.stats.multitest import multipletests



def load_data_from_csv(infile_path):
    header=np.genfromtxt(infile_path, delimiter=",", max_rows=1, dtype='str')
    data=np.genfromtxt(infile_path, delimiter=",", skip_header=1, dtype='str')
    gene_id = data[:,0]
    headers = list(header)
    arr=data[:,1:]
    arr_clean = np.where(np.isin(arr, ["NA", "N/A", "", "--"]), "nan", arr)
    numerical_data = arr_clean.astype(float)
    #print(normalized_data.shape)
    print(f"Successfully loaded {infile_path}")
    return headers, gene_id, numerical_data