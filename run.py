#! /data/keeling/a/nmd/miniconda3/envs/sage_full/bin/sage-python -u

#SBATCH --array=0-19
#SBATCH --partition m
#SBATCH --tasks=1
#SBATCH --mem-per-cpu=4G
#SBATCH --nice=10000
#SBATCH --time=7-00:00
#SBATCH --output=/data/keeling/a/chaeryn2/SLURM_print/count_components%A_%a
#SBATCH --error=/data/keeling/a/chaeryn2/SLURM_error/count_componenets%A_%a

from count_components_manifold import *
from orbits_manifold import *
import os
import multiprocessing
import gc
import pandas as pd
import snappy.snap.t3mlite as t3m

def ahg_for_manifolds():
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

def ahg_randomize(M):
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
        p = multiprocessing.Process(target=ahg_randomize(), args=M)
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