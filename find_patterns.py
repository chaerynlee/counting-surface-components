from count_components_manifold import *
from orbits_manifold import *
from Orbits import orbits
import os
import pandas as pd
import snappy.snap.t3mlite as t3m
import itertools
import random
import math


# TODO
# - modify simplify_remove_one() to check more cases (not just removing the first insignificant pairing in the list but maybe all?)
# - write new function that checks all combinations of pairings of a certain small number(maybe 5?)
# - run both functions on an example to see if we get reasonable results analyzable by hand

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

if __name__ == '__main__':
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
