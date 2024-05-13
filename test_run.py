from count_components_manifold import *
from orbits_manifold import *
from Orbits import orbits
import os
import pandas as pd
import snappy.snap.t3mlite as t3m

def aht_for_manifolds():
    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    index = df.index[df['name'] == 't12199'].values[0]
    pattern = df['gen_func'].iloc[index]
    df_pattern = df[df['gen_func'] == pattern]
    mflds = df_pattern['name'].tolist()

    # ran until index 166 (needed higher euler bound for ConnectedSurfaces)
    # was successful for 'o9_44238'
    for name in mflds:
        print(name)
        M = snappy.Manifold(name)
        CS = ConnectedSurfaces(M)
        LW = CS.essential_faces_of_normal_polytope()
        LW_faces = LW.maximal
        for i in range(len(LW_faces)):
            vs = LW_faces[i].vertex_surfaces
            vs_regina_list = [S.surface for S in vs]
            SO = SurfacetoOrbit(vs_regina_list)
            G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
            G.reduce()

    # with open('example_default_tri.pickle', 'rb') as f:
    #     eg_int, ed_int_div, eg_pairings = pickle.load(f)
    #
    # G = Pseudogroup(eg_pairings, eg_int, ed_int_div)
    # G.reduce()

def aht_randomize(M):
    # M: snappy manifold
    tri_found = False
    tri_isosig = []
    while not tri_found:
        CS = ConnectedSurfaces(M, -10)
        LW = CS.essential_faces_of_normal_polytope()
        LW_faces = LW.maximal

        for i in range(len(LW_faces)):
            vs = LW_faces[i].vertex_surfaces
            vs_regina_list = [S.surface for S in vs]
            SO = SurfacetoOrbit(vs_regina_list)
            G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
            result = G.reduce()
            if result:
                return result
        tri_isosig.append(M.triangulation_isosig())
        print('Not found')
        while M.triangulation_isosig() in tri_isosig:
            M.randomize()
            print('Randomizing')

def print_manifold_info(M):
    # M: snappy manifold
    print('num tet:', M.num_tetrahedra())
    Mcomplex = t3m.Mcomplex(M)
    print()
    CS = ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    LW_faces = LW.maximal
    for i in range(len(LW_faces)):
        vs = LW_faces[i].vertex_surfaces
        vs_regina_list = [S.surface for S in vs]
        for S in vs_regina_list:
            print(S.detail())

def example():
    M = snappy.Manifold('K13n586')
    CS = ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    LW_faces = LW.maximal
    vs = LW_faces[0].vertex_surfaces
    vs_regina_list = [S.surface for S in vs]
    SO = SurfacetoOrbit(vs_regina_list)
    return SO.pairings[0]


if __name__ == '__main__':
    # Manifolds that are known to have the given pattern
    # -Monetesinos knots with exactly 4 fractions (rational tangle) and where all denominators are larger than or equal to 3

    # -In KnotInfo found 12a_554, 12a_740, 12n_553, 12n_554, 12n_555, 12n_556 (should filter for more)
    # -Check Nathan's paper for complicated examples


    M = snappy.Manifold('K13n586_nice.tri')
    CS = ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    LW_faces = LW.maximal
    for i in range(len(LW_faces)):
        vs = LW_faces[i].vertex_surfaces
        vs_regina_list = [S.surface for S in vs]
        SO = SurfacetoOrbit(vs_regina_list)

        assign = {'x0':2, 'x1':3}

        assigned_interval = orbits.Interval(SO.interval.start.evaluate(**assign) + 1, SO.interval.end.evaluate(**assign))
        print('interval', assigned_interval)
        assigned_pairings = []
        for p in SO.pairings:
            domain_start = p.domain.start.evaluate(**assign)
            domain_end = p.domain.end.evaluate(**assign)
            iso_shift = p.isometry.shift.evaluate(**assign)
            iso_flip = p.isometry.flip
            if not iso_flip:
                assigned_pairings.append(orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift, iso_flip)))
            if iso_flip:
                assigned_pairings.append(orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift + 1, iso_flip)))
        print('pairings')
        for p in sorted(assigned_pairings):
            print(p)

        G = orbits.Pseudogroup(assigned_pairings, assigned_interval)
        print(G.reduce())

        # G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
        # result = G.reduce()
        # print(result)
