from count_components_manifold import *
from orbits_manifold import *
from nscomplex_updated import *
import count_components
import os, itertools, random, math, time, copy
import pandas as pd

def significant_col(interval, pairings, num_var, subcollection_indices, var_range=50):
    """
    Returns True or False depending on whether the given subcollection of pairings is significant.
    By default, evaluates the original collection and subcollection of pairings for all variables up to 50.
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

# def extend_gen_fcn(manifold, genus):
#     """
#     For the given manifold, calculates the number of connected surfaces of the given genus.
#     manifold must be a string of its name and must be contained in the data set 'very_large_combined.csv'
#     """
#     df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
#     i = df.index[df['name'] == manifold].values[0]
#     TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
#     T = regina.Triangulation3(TS)
#     LWC_info = df.iloc[i, df.columns.get_loc('all_faces')]
#     vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
#     genera_info = df.iloc[i, df.columns.get_loc('vertex_genera')]
#
#     count = 0
#     for face in eval(LWC_info):
#         # print('face', face)
#         surface_names = face['verts']
#         dim = face['dim']
#
#         vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
#         # print('vertex surfaces', vertex_surface_vectors)
#         vertex_surface_genera = [eval(genera_info)[name] for name in surface_names]
#         vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
#         vertex_surfaces_ns = [surfaces.NormalSurface(S, i) for i, S in enumerate(vertex_surfaces)]
#         num_var = len(surface_names)
#
#         AF = faces.AdmissibleFace(dim, vertex_surfaces_ns)
#         solutions = AF.surfaces_of_potential_genus_in_interior(genus)
#         # print(len(solutions))
#
#         for sol in solutions:
#             # print(sol)
#             # print(sol.isConnected())
#             num_comp = count_components.SurfacetoOrbit(sol).countcomponents()
#             # print(num_comp)
#             if num_comp == 1:
#                 count += 1
#     return count

def extend_gen_fcn(manifold, genus, all=False):
    """
    For the given manifold, calculates the number of connected surfaces of the given genus.
    manifold must be a string of its name and must be contained in the data set 'very_large_combined.csv'
    If all is set to True returns a list of all numbers of connected surfaces up to the given genus.
    Also returns the time taken to perform the computation.
    """
    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == manifold].values[0]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    LWC_info = df.iloc[i, df.columns.get_loc('all_faces')]
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    genera_info = df.iloc[i, df.columns.get_loc('vertex_genera')]

    tik = time.perf_counter()
    if not all:
        count = count_conn_surfaces(LWC_info, T, genera_info, genus, vector_info)
        tok = time.perf_counter()
        time_taken = tok - tik
        return count, time_taken
    if all:
        counts = []
        for g in range(2, genus + 1):
            counts.append(count_conn_surfaces(LWC_info, T, genera_info, g, vector_info))
        tok = time.perf_counter()
        time_taken = tok - tik
        return counts, time_taken

def count_conn_surfaces(LWC_info, T, genera_info, genus, vector_info):
    """
    Helper function for extend_gen_fcn that constructs vertex surfaces from the given information for each LW-face,
    finds all potential surfaces of the given genus then determines whether it is indeed connected.
    Returns the count of all connected surfaces of the given genus.
    """
    count = 0
    for face in eval(LWC_info):
        # print('face', face)
        surface_names = face['verts']
        dim = face['dim']

        vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
        # print('vertex surfaces', vertex_surface_vectors)
        vertex_surface_genera = [eval(genera_info)[name] for name in surface_names]
        vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
        vertex_surfaces_ns = [surfaces.NormalSurface(S, i) for i, S in enumerate(vertex_surfaces)]
        num_var = len(surface_names)

        AF = faces.AdmissibleFace_nozeroset(dim, vertex_surfaces_ns)
        solutions = AF.surfaces_of_potential_genus_in_interior(genus)
        # print(len(solutions))

        for sol in solutions:
            # print(sol)
            # print(sol.isConnected())
            num_comp = count_components.SurfacetoOrbit(sol).countcomponents()
            # print(num_comp)
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
    # test all manifolds to check if extend_gen_fcn is correct
    # df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    # for i in range(1, len(df.index)):
    #     M = df.iloc[i, df.columns.get_loc('name')]
    #     actual_count = eval(df.iloc[i, df.columns.get_loc('by_genus')])
    #     count = extend_gen_fcn(M, 21, all=True)
    #     print(M, actual_count == count)

    M = 'K15n129923'

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == M].values[0]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]

    transforms = [[{'t0': 1, 't1': 0}, {'t0': 3, 't1': 1}],
                  [{'t0': 1, 't1': 1}, {'t0': 2, 't1': 3}],
                  [{'t0': 2, 't1': 1}, {'t0': 3, 't1': 2}],
                  [{'t0': 1, 't1': 2}, {'t0': 1, 't1': 3}],
                  [{'t0': 3, 't1': 1}, {'t0': 2, 't1': 1}],
                  [{'t0': 2, 't1': 3}, {'t0': 1, 't1': 2}],
                  [{'t0': 3, 't1': 2}, {'t0': 1, 't1': 1}],
                  [{'t0': 1, 't1': 3}, {'t0': 0, 't1': 1}]]

    # transforms = [[{'t0': 1, 't1': 0}, {'t0': 1, 't1': 1}],
    #               [{'t0': 1, 't1': 1}, {'t0': 0, 't1': 1}]]

    # inverse_transforms = [[{'x0': 1, 'x1': 0}, {'x0': -3, 'x1': 1}],
    #                       [{'x0': 3, 'x1': -1}, {'x0': -2, 'x1': 1}],
    #                       [{'x0': 2, 'x1': -1}, {'x0': -3, 'x1': 2}],
    #                       [{'x0': 3, 'x1': -2}, {'x0': -1, 'x1': 1}],
    #                       [{'x0': 1, 'x1': -1}, {'x0': -2, 'x1': 3}],
    #                       [{'x0': 2, 'x1': -3}, {'x0': -1, 'x1': 2}],
    #                       [{'x0': 1, 'x1': -2}, {'x0': -1, 'x1': 3}],
    #                       [{'x0': 1, 'x1': -3}, {'x0': 0, 'x1': 1}]]

    # inverse_transforms = [[{'x0': 1, 'x1': 0}, {'x0': 1, 'x1': -1}],
    #                       [{'x0': 1, 'x1': -1}, {'x0': 0, 'x1': 1}]]

    for face_num, face in enumerate(eval(LWC_info)):
        print('face', face_num)
        surface_names = face['verts']
        vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
        vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in
                           vertex_surface_vectors]
        # SO = SurfacetoOrbit(vertex_surfaces)
        # X = vertex_surfaces[0]
        # Y = vertex_surfaces[1]
        # G = Pseudogroup_comparable(SO.pairings, SO.interval)
        #
        # for comb in itertools.product(range(1, 10), repeat=2):
        #     assign = dict()
        #     for n in range(2):
        #         var_name = 'x' + str(n)
        #         assign[var_name] = comb[n]
        #     # S = count_components.SurfacetoOrbit(X*regina.LargeInteger(assign['x0']) + Y*regina.LargeInteger(assign['x1']))
        #     G_num = evaluate_pseudogroup(G.universe, G.pairings, assign)
        #     print(assign, G_num)

        for j, t in enumerate(transforms):
            print('case', j, t)
            SO = SurfacetoOrbit(vertex_surfaces)
            G = Pseudogroup_comparable(SO.pairings, SO.interval)
            G.transform(t)

            count = G.reduce()
            print(count)
            # for comb in itertools.product(range(1, 10), repeat=2):
            #     assign = dict()
            #     for n in range(2):
            #         var_name = 't' + str(n)
            #         assign[var_name] = comb[n]
            #     G_num = evaluate_pseudogroup(G.universe, G.pairings, assign)
            #     print(assign, G_num)
            # if (isinstance(count, Polynomial)):
            #     print(linear_transform_poly(count, inverse_transforms[j]))
            print()


    # for n in range(1, 7):
    #     for m in range(1, 7):
    #         count = evaluate_pseudogroup(SO.interval, SO.pairings, {'x0': n, 'x1': m})
    #         S = vertex_surfaces[0] * regina.LargeInteger(n) + vertex_surfaces[1] * regina.LargeInteger(m)
    #         print(n, m, S.eulerChar(), S.isConnected(), count_components.SurfacetoOrbit(S).countcomponents(), count)
    # print()

    # M = snappy.Manifold('jLLPzPQcdeffghiiihsteaviivg_bBBa')
    # print(M.triangulation_isosig())
    # CS = connected_surfaces.ConnectedSurfaces(M)
    # LW = CS.essential_faces_of_normal_polytope()
    # LW_max = LW.maximal
    # print(len(LW_max))
    # for face in LW_max:
    #     print(face)
    #     vertex_surfaces = [S.surface for S in face.vertex_surfaces]
    #     SO = SurfacetoOrbit(vertex_surfaces)
    #     G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
    #     simplified_interval, simplified_pairings = G.reduce_amap()
    #
    #     found = False
    #     for n in range(2, 4):
    #         result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
    #         if result:
    #             print('simplified_interval', simplified_interval)
    #             print('pattern', result)
    #             found = True
    #             break
    #
    #     if not found:
    #         print('not found for some face')
    #         break
    #     else:
    #         continue

