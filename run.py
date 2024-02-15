from count_components_manifold import *
from orbits_manifold import *
import os
import pandas as pd

if __name__ == '__main__':
    df = pd.read_csv(os.getcwd() + '\cusped_census_very_large.csv')
    index = df.index[df['name'] == 't12198']
    pattern = df['gen_func'].iloc[31]
    df_pattern = df[df['gen_func'] == pattern]
    mflds = df_pattern['name'].tolist()

    with open('example_default_tri.pickle', 'rb') as f:
        eg_int, ed_int_div, eg_pairings = pickle.load(f)

    G = Pseudogroup(eg_pairings, eg_int, ed_int_div)
    G.reduce()