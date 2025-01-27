#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=0-49
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/count_components%A_%a
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/count_components%A_%a

from count_components_manifold import *
from orbits_manifold import *
from find_patterns import *
import os
import multiprocessing
import gc
import pandas as pd
import snappy.snap.t3mlite as t3m
import copy

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
                vertex_surfaces = [[S.full_vector for S in LW_faces[i].vertex_surfaces] for i in range(len(LW_faces))]
                save = {'manifold': M.name(),
                        'triangulation': M.triangulation_isosig(),
                        'LW_complex': vertex_surfaces,
                        'candidates': result}
                directory = '/data/keeling/a/chaeryn2/ahg_results/'
                filename = f'PG_info_{M.name()}'
                with open(directory+filename, 'wb') as file:
                    pickle.dump(save, file)
                return
        tri_isosig.append(M.triangulation_isosig())
        while M.triangulation_isosig() in tri_isosig:
            M.randomize()

def find_pattern_unknown(index=-1):
    df_all = pd.read_csv(os.getcwd() + '/very_large_combined.csv')

    f = open(os.getcwd() + '/manifolds_with_least_LWfaces.txt')
    mflds = f.read().split('\n')
    df = df_all[df_all['name'].isin(mflds)]

    tri_info = df['tri_used'].tolist()
    vector_info = df['vertex_surfaces'].tolist()
    LWC_info = df['max_faces'].tolist()

    if index == -1:
        for i, M in df['name'].tolist():
            TS = snappy.Manifold(tri_info[i])
            T = regina.Triangulation3(TS)
            interval_allfaces = []
            result_allfaces = []
            original_PG = []
            for face in eval(LWC_info[i]):
                surface_names = face['verts']
                if len(surface_names) == 1:
                    interval_allfaces.append('single_vertex_surface')
                    result_allfaces.append('single_vertex_surface')
                else:
                    vertex_surface_vectors = [eval(vector_info[i])[name] for name in surface_names]
                    vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in
                                       vertex_surface_vectors]
                    SO = SurfacetoOrbit(vertex_surfaces)
                    G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
                    original_PG.append(G)
                    simplified_interval, simplified_pairings = G.reduce_amap()
                    interval_allfaces.append(simplified_interval)

                    # test all subcollections of size 2-6, stop if something is found
                    n = 2
                    for n in range(2, 7):
                        result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
                        if result:
                            break
                        else:
                            continue
                    # for time efficiency just return 'not_found' if the previous step fails
                    if not result:
                        # result = simplify_remove_one(simplified_interval, simplified_pairings, SO.num_vertex)
                        result = 'not_found'
                    result_allfaces.append(result)

            save = {'manifold': M,
                    'LW_complex': LWC_info[i],
                    'orginal_psuedogroup': original_PG,
                    'intervals': interval_allfaces,
                    'patterns': result_allfaces}
            filename = f'unknown_pattern_info_{M}'
            with open(filename, 'wb') as file:
                pickle.dump(save, file)
    else:
        i = index
        M = mflds[i]
        TS = snappy.Manifold(tri_info[i])
        T = regina.Triangulation3(TS)
        interval_allfaces = []
        result_allfaces = []
        original_PG = []
        for face in eval(LWC_info[i]):
            surface_names = face['verts']
            if len(surface_names) == 1:
                interval_allfaces.append('single_vertex_surface')
                result_allfaces.append('single_vertex_surface')
            else:
                vertex_surface_vectors = [eval(vector_info[i])[name] for name in surface_names]
                vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in
                                   vertex_surface_vectors]
                SO = SurfacetoOrbit(vertex_surfaces)
                G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
                original_PG.append(G)
                simplified_interval, simplified_pairings = G.reduce_amap()
                interval_allfaces.append(simplified_interval)

                # test all subcollections of size 2-6, stop if something is found
                n = 2
                for n in range(2, 7):
                    result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
                    if result:
                        break
                    else:
                        continue
                # for time efficiency just return 'not_found' if the previous step fails
                if not result:
                    # result = simplify_remove_one(simplified_interval, simplified_pairings, SO.num_vertex)
                    result = 'not_found'
                result_allfaces.append(result)

        save = {'manifold': M,
                'LW_complex': LWC_info[i],
                'orginal_psuedogroup': original_PG,
                'intervals': interval_allfaces,
                'patterns': result_allfaces}
        filename = f'unknown_pattern_info_{M}'
        with open(filename, 'wb') as file:
            pickle.dump(save, file)


def find_pattern(M):
    correct_euler = False
    euler_bound = -6
    while not correct_euler:
        try:
            CS = ConnectedSurfaces(M, euler_bound)
            correct_euler = True
            LW = CS.essential_faces_of_normal_polytope()
            LW_faces = LW.maximal
        except:
            correct_euler = False
            euler_bound += -2

    result_allfaces = []
    interval_allfaces = []
    for i in range(len(LW_faces)):
        vs = LW_faces[i].vertex_surfaces
        vs_regina_list = [S.surface for S in vs]
        SO = SurfacetoOrbit(vs_regina_list)
        G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
        simplified_interval, simplified_pairings = G.reduce_amap()
        interval_allfaces.append(simplified_interval)

        # test all subcollections of size 2-6, stop if something is found
        n = 2
        for n in range(2, 7):
            result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
            if result:
                break
            else:
                continue
        # if no significant subcollection of size at most 6 is not found, simplify by removing pairings one at a time
        if not result:
            result = simplify_remove_one(simplified_interval, simplified_pairings, SO.num_vertex)
        result_allfaces.append(result)

    vertex_surfaces = [[S.full_vector for S in LW_faces[i].vertex_surfaces] for i in range(len(LW_faces))]
    save = {'manifold': M.name(),
            'LW_complex': vertex_surfaces,
            'intervals': interval_allfaces,
            'patterns': result_allfaces}
    directory = '/data/keeling/a/chaeryn2/patterns/'
    filename = f'pattern_info_{M.name()}'
    with open(directory + filename, 'wb') as file:
        pickle.dump(save, file)
    return

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

def find_smallest_mfld_by_genfcn():
    # this function has already been run and the necessary txt file has been made, just left in case of future use
    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    gen_fcn = df['gen_func'].unique().tolist()
    least_time = []
    for gf in gen_fcn:
        df_gf = df[df['gen_func'] == gf]
        index = df_gf['gen_func_time'].idxmin()
        least_time.append(df_gf.loc[index, 'name'])

    with open('smallest_manifold_by_genfcn.txt', 'w') as f:
        for name in least_time:
            f.write(str(name) + '\n')

def find_rep_manifold_ebg():
    df = pd.read_csv(os.getcwd() + '/extended_by_genus.csv')
    df_all = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    gen_fcn = df['gen_func'].unique().tolist()
    least_faces = []
    for gf in gen_fcn:
        df_gf = df_all[df_all['gen_func'] == gf]
        min_faces = df_gf['LW_num_max_faces'].min()
        df_gf_mf = df_gf[df_gf['LW_num_max_faces'] == min_faces]
        index = df_gf_mf['tets'].idxmin()
        least_faces.append(df_gf_mf.loc[index, 'name'])

    with open('manifolds_with_least_LWfaces.txt', 'w') as f:
        for name in least_faces:
            f.write(str(name) + '\n')

def main_aht_randomize():
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    index = df.index[df['name'] == 't12199'].values[0]
    pattern = df['gen_func'].iloc[index]
    df_pattern = df[df['gen_func'] == pattern]
    mflds = df_pattern['name'].tolist()

    mfld_list = []
    with open(f'exceeded_time.txt', 'r') as f:
        exceeded_time = f.read().split()

    for i in range(task, len(mflds), 20):
        found = False
        name = mflds[i]
        for filename in os.listdir('/data/keeling/a/chaeryn2/ahg_results/'):
            if name in filename:
                found = True
                break
        if not found:
            if i not in exceeded_time:
                mfld_list.append(name)

    for name in mfld_list:
        gc.collect()
        M = snappy.Manifold(name)
        p = multiprocessing.Process(target=aht_randomize, args=(M,))
        p.start()
        p.join(5000)
        if p.is_alive():
            flag = True
            while flag:
                if 'dummy_file' not in os.listdir():
                    with open('dummy_file', 'w') as f:
                        pass
                    with open(f'exceeded_time.txt', 'a') as file:
                        file.write("\n" + M.name() + "\n")
                    os.remove('dummy_file')
                    p.terminate()
                    flag = False
                else:
                    continue
            continue

def main_find_pattern():
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    mflds = df['name'].tolist()

    mfld_list = []
    for i in range(task, len(mflds), 20):
        found = False
        name = mflds[i]
        for filename in os.listdir('/data/keeling/a/chaeryn2/patterns/'):
            if name in filename:
                found = True
                break
        if not found:
            mfld_list.append(name)

    for name in mfld_list:
        M = snappy.Manifold(name)
        find_pattern(M)

def main_find_pattern_unknown():
    i = int(os.environ['SLURM_ARRAY_TASK_ID'])

    df_all = pd.read_csv(os.getcwd() + '/very_large_combined.csv')

    f = open(os.getcwd() + '/manifolds_with_least_LWfaces.txt')
    mflds = f.read().split('\n')
    df = df_all[df_all['name'].isin(mflds)]

    # mfld_list = []
    # for i in range(task, df.shape[0], 76):
    #     found = False
    #     name = mflds[i]
    #     for filename in os.listdir('/data/keeling/a/chaeryn2/patterns/'):
    #         if name in filename:
    #             found = True
    #             break
    #     if not found:
    #         mfld_list.append(name)

    M = df.iloc[i, df.columns.get_loc('name')]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]

    interval_allfaces = []
    result_allfaces = []
    original_PG = []
    for face in eval(LWC_info):
        surface_names = face['verts']
        if len(surface_names) == 1:
            interval_allfaces.append('single_vertex_surface')
            result_allfaces.append('single_vertex_surface')
        else:
            vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
            vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
            SO = SurfacetoOrbit(vertex_surfaces)
            G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
            G_copy = copy.deepcopy(G)
            original_PG.append(G_copy)
            simplified_interval, simplified_pairings = G.reduce_amap()
            interval_allfaces.append(simplified_interval)

            # test all subcollections of size 2-6, stop if something is found
            for n in range(2, 7):
                result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
                if result:
                    break
                else:
                    continue
            # for time efficiency just return 'not_found' if the previous step fails
            if not result:
                # result = simplify_remove_one(simplified_interval, simplified_pairings, SO.num_vertex)
                result = 'not_found'
            result_allfaces.append(result)

    save = {'manifold': M,
            'LW_complex': LWC_info,
            'orginal_psuedogroup': original_PG,
            'intervals': interval_allfaces,
            'patterns': result_allfaces}
    directory = '/data/keeling/a/chaeryn2/patterns/'
    filename = f'unknown_pattern_info_{M}'
    with open(directory + filename, 'wb') as file:
        pickle.dump(save, file)

def main_find_pattern_unknown_separate():
    i = int(os.environ['SLURM_ARRAY_TASK_ID'])

    df_all = pd.read_csv(os.getcwd() + '/very_large_combined.csv')

    f = open(os.getcwd() + '/manifolds_with_least_LWfaces.txt')
    mflds = f.read().split('\n')
    df = df_all[df_all['name'].isin(mflds)]
    M = df.iloc[i, df.columns.get_loc('name')]

    # check if manifold already has results
    for filename in os.listdir('/data/keeling/a/chaeryn2/patterns/'):
        if M in filename:
            return

    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]

    for face_num, face in enumerate(eval(LWC_info)):
        surface_names = face['verts']
        if len(surface_names) == 1:
            continue
        else:
            vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
            vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
            SO = SurfacetoOrbit(vertex_surfaces)
            G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
            G_copy = copy.deepcopy(G)

            simplified_interval, simplified_pairings = G.reduce_amap()

            # test all subcollections of size 2-6, stop if something is found
            for n in range(2, 7):
                result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
                if result:
                    save = {'manifold': M,
                            'LW_complex': LWC_info,
                            'orginal_psuedogroup': G_copy,
                            'intervals': simplified_interval,
                            'patterns': result}
                    directory = '/data/keeling/a/chaeryn2/patterns/'
                    filename = f'unknown_pattern_info_{M}_face{face_num}_subcoll{n}'
                    with open(directory + filename, 'wb') as file:
                        pickle.dump(save, file)
                    break
                else:
                    save = {'manifold': M,
                            'LW_complex': LWC_info,
                            'orginal_psuedogroup': G_copy,
                            'intervals': simplified_interval,
                            'patterns': 'not_found'}
                    directory = '/data/keeling/a/chaeryn2/patterns/'
                    filename = f'unknown_pattern_info_{M}_face{face_num}_subcoll{n}'
                    with open(directory + filename, 'wb') as file:
                        pickle.dump(save, file)

def main_find_pattern_by_genfcn():
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    with open('smallest_manifold_by_genfcn.txt', 'r') as f:
        mflds = f.read().split('\n')

    mfld_list = []
    for i in range(task, len(mflds), 20):
        found = False
        name = mflds[i]
        for filename in os.listdir('/data/keeling/a/chaeryn2/patterns/'):
            if name in filename:
                found = True
                break
        if not found:
            mfld_list.append(name)

    for name in mfld_list:
        M = snappy.Manifold(name)
        find_pattern(M)

def recreate_example(M):
    # M: string of the manifold's name
    print('manifold', M)

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == M].values[0]
    # M = df.iloc[i, df.columns.get_loc('name')]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]

    interval_allfaces = []
    result_allfaces = []
    original_PG = []
    for face in eval(LWC_info):
        print(face)
        surface_names = face['verts']
        if len(surface_names) == 1:
            interval_allfaces.append('single_vertex_surface')
            result_allfaces.append('single_vertex_surface')
        else:
            vertex_surface_vectors = [eval(vector_info)[name] for name in surface_names]
            vertex_surfaces = [regina.NormalSurface(T, regina.NS_QUAD_CLOSED, vec) for vec in vertex_surface_vectors]
            SO = SurfacetoOrbit(vertex_surfaces)
            G = Pseudogroup(SO.pairings, SO.interval, SO.interval_divided)
            print(G)
            original_PG.append(G)
            simplified_interval, simplified_pairings = G.reduce_amap()
            interval_allfaces.append(simplified_interval)
            print('simplified to', simplified_interval)
            print(simplified_pairings)

            # test all subcollections of size 2-6, stop if something is found
            for n in range(2, 7):
                result = test_all_subcol(simplified_interval, simplified_pairings, SO.num_vertex, n)
                if result:
                    break
                else:
                    continue
            # for time efficiency just return 'not_found' if the previous step fails
            if not result:
                # result = simplify_remove_one(simplified_interval, simplified_pairings, SO.num_vertex)
                result = 'not_found'
            print('result')
            for r in result:
                print(r)
            result_allfaces.append(result)

    save = {'manifold': M,
            'LW_complex': LWC_info,
            'orginal_psuedogroup': original_PG,
            'intervals': interval_allfaces,
            'patterns': result_allfaces}
    filename = f'unknown_pattern_info_{M}'
    with open(filename, 'wb') as file:
        pickle.dump(save, file)

def check_extend_gen_fcn():
    task = int(os.environ['SLURM_ARRAY_TASK_ID'])

    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')

    for i in range(910+task, len(df.index), 50):
        M = df.iloc[i, df.columns.get_loc('name')]
        actual_count = eval(df.iloc[i, df.columns.get_loc('by_genus')])
        count = extend_gen_fcn(M, 21, all=True)
        save = {'manifold': M,
                'count': count,
                'actual_count': actual_count,
                'match': actual_count == count}

        directory = '/data/keeling/a/chaeryn2/patterns/'
        filename = f'gen_fcn_{M}'
        with open(directory + filename, 'wb') as file:
            pickle.dump(save, file)

if __name__ == '__main__':
    check_extend_gen_fcn()
    # main_find_pattern_unknown_separate()
    # recreate_example('o9_34491')