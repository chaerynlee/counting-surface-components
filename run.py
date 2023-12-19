from count_components_manifold import *
from orbits_manifold import *

if __name__ == '__main__':
    # unknown is considered as True
    with open('example.pickle', 'rb') as f:
        eg_int, eg_pairings = pickle.load(f)

    G = Pseudogroup(eg_pairings, eg_int)
    G.reduce()