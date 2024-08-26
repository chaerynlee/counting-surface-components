from count_components_manifold import *
from orbits_manifold import *
from Orbits import orbits
import count_components
import os, itertools, random, math
import pandas as pd
import numpy as np
import snappy.snap.t3mlite as t3m

def evaluate_interval(interval, assign):
    return orbits.Interval(interval.start.evaluate(**assign) + 1, interval.end.evaluate(**assign))

def evaluate_pairing(pairing, assign):
    domain_start = pairing.domain.start.evaluate(**assign)
    domain_end = pairing.domain.end.evaluate(**assign)
    iso_shift = pairing.isometry.shift.evaluate(**assign)
    iso_flip = pairing.isometry.flip
    if not iso_flip:
        return orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift, iso_flip))
    if iso_flip:
        return orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift + 1, iso_flip))

def evaluate_pseudogroup(interval, pairings, assign):
    """
    Evaluates given interval and pairings of polynomials with the variables in assign (must be a dictionary with the right variable names)
    and returns the number of orbits.
    """
    assigned_interval = orbits.Interval(interval.start.evaluate(**assign) + 1, interval.end.evaluate(**assign))
    assigned_pairings = []
    for p in pairings:
        domain_start = p.domain.start.evaluate(**assign)
        domain_end = p.domain.end.evaluate(**assign)
        iso_shift = p.isometry.shift.evaluate(**assign)
        iso_flip = p.isometry.flip
        if not iso_flip:
            assigned_pairings.append(orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift, iso_flip)))
        if iso_flip:
            assigned_pairings.append(orbits.Pairing(orbits.Interval(domain_start + 1, domain_end), orbits.Isometry(iso_shift + 1, iso_flip)))

    G = orbits.Pseudogroup(assigned_pairings, assigned_interval)
    return G.reduce()

def significant_col(interval, pairings, num_var, subcollection_indices, var_range=50):
    """
    Returns True or False depending on whether the given subcollection of pairings is significant.
    The subcollection is given as their indices not as pairings themselves.
    """
    subcollection = [pairings[i] for i in subcollection_indices]
    for comb in itertools.product(range(1, var_range), repeat=num_var):
        assign = dict()
        for n in range(num_var):
            var_name = 'x' + str(n)
            assign[var_name] = comb[n]
        actual = evaluate_pseudogroup(interval, pairings, assign)
        sub = evaluate_pseudogroup(interval, subcollection, assign)

        if actual != sub:
            # print(assign)
            # print('actual num components', actual)
            # print('num components of subcollection', sub)
            return False
        else:
            continue
    return True

def simplify_halving(interval, pairings, num_var):
    """
    Attempts to simplify given pseudogroup by taking half of the pairings.
    This algorithm will continue to take the half of whatever subcollection of pairings it is successful in checking.
    Randomizes the order of pairings at each attempt.
    """
    random.shuffle(pairings)

    masterlist = [list(range(len(SO.pairings)))]
    masterlist_before = []
    # [SO.pairings] -> [A, B] -> [A1, B1, B2] -> [A11, A12, B21] -> []
    while masterlist != []:
        update_masterlist = []
        for subcol in masterlist:
            midpoint = math.ceil(len(subcol) / 2)
            subcol1 = subcol[0:midpoint]
            subcol2 = subcol[midpoint:]
            result1 = significant_col(interval, pairings, num_var, subcol1)
            result2 = significant_col(interval, pairings, num_var, subcol2)
            if result1:
                update_masterlist.append(subcol1)
            if result2:
                update_masterlist.append(subcol2)
        masterlist_before = masterlist
        masterlist = update_masterlist

    if masterlist_before == [list(range(len(SO.pairings)))]:
        return False
    else:
        sub_pairings = []
        for sublist in masterlist_before:
            sub_pairings.append([pairings[i] for i in sublist])
        return sub_pairings

def simplify_remove_one(interval, pairings, num_var):
    """
    Remove the first insignificant pairing in the list of pairings.
    Do this as many times as necessary until it gets to the smallest set of pairings possible.
    """
    original_len = len(pairings)
    if original_len < 4:
        return pairings

    simplify_possible = True
    while simplify_possible:
        simplify_possible = False
        for i, p in enumerate(pairings):
            # print(i, p)
            l = list(range(len(pairings)))
            l.remove(i)
            result = significant_col(interval, pairings, num_var, l)
            # print(result)
            if result:
                pairings.remove(p)
                simplify_possible = True
                break
    if len(pairings) == original_len:
        return False
    else:
        return pairings

def test_all_subcol(interval, pairings, num_var, n=5):
    """
    Takes all subcollections of size n and returns the first subcollection that is possibly significant.
    By default, checks subcollections of size 5.
    """
    for comb in itertools.combinations(range(len(pairings)), n):
        result = significant_col(interval, pairings, num_var, comb)
        if result:
            subcollection = [pairings[i] for i in comb]
            return subcollection
    return False

def extend_gen_fcn(manifold, genus):
    """
    For the given manifold, calculates the number of connected surfaces of the given genus.
    manifold must be a string of its name and must be contained in the data set 'very_large_combined.csv'
    """
    # print('manifold', manifold)

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == manifold].values[0]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    genera_info = df.iloc[i, df.columns.get_loc('vertex_genera')]

    count = 0
    for face in eval(LWC_info):
        # print('face', face)
        surface_names = face['verts']
        if len(surface_names) == 1:
            # print('single vertex surface')
            surface_genus = eval(genera_info)[surface_names[0]]
            if surface_genus == genus:
              count += 1
        else:
            # print('more than one vertex surface')
            vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
            # print('vertex surfaces', vertex_surface_vectors)
            vertex_surface_genera = [eval(genera_info)[name] for name in surface_names]
            vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
            num_var = len(surface_names)

            # TODO
            # is there a quick way to find all surfaces of given genus (all lattice points in cone of each lw-face)
            # from just the vertex surface vectors?
            solutions = [(1, 2), (3, 4), (5, 6)]  # just some random example to test code, should be replaced once TODO is solved

            for sol in solutions:
                sol_vector = np.zeros(len(vertex_surface_vectors[0]))
                for i, k in enumerate(sol):
                    sol_vector += k * np.array(vertex_surface_vectors[i])
                regina_vec = [regina.LargeInteger(n) for n in sol_vector]
                sol_surface = regina.NormalSurface(T, regina.NS_QUAD_CLOSED, regina_vec)
                num_comp = count_components.SurfacetoOrbit(sol_surface).countcomponents()
                if num_comp == 1:
                    count += 1
    return count


def some_tests():
    M = snappy.Manifold('K13n586_nice.tri')
    CS = ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    LW_faces = LW.maximal
    for i in range(len(LW_faces)):
        vs = LW_faces[i].vertex_surfaces
        vs_regina_list = [S.surface for S in vs]
        SO = SurfacetoOrbit(vs_regina_list)
        # print('unsimplified')
        # print('interval', SO.interval)
        # print('pairings', len(SO.pairings))
        # for p in SO.pairings:
        #     print(p)
        # print()
        #
        # result_ver1 = simplify_remove_one(SO.interval, SO.pairings, 2)
        # print('removing one', result_ver1)
        #
        # print()
        #
        # result_ver2 = test_all_subcol(SO.interval, SO.pairings, 2)
        # print('testing subcollections', result_ver2)

        print('simplified')
        G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
        simplified_interval, simplified_pairings = G.reduce_amap()

        print('interval', simplified_interval)
        print('pairings', len(simplified_pairings))
        for p in simplified_pairings:
            print(p)
        print()

        # result_ver1 = simplify_remove_one(simplified_interval, simplified_pairings, 2)
        # print('removing one', result_ver1)

        print()

        result_ver2 = test_all_subcol(simplified_interval, simplified_pairings, 2, 6)
        print('testing subcollections', result_ver2)

        # assign = {'x0': 2, 'x1': 3}
        # result = simplify_remove_one(SO)
        # print('result', result)
        # result = test_all_subcol(SO, 3)
        # print(result)


if __name__ == '__main__':
    for n in range(2, 22):
        print(extend_gen_fcn('s783', n))
