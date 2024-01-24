import snappy, regina, nscomplex
import snappy.snap.t3mlite as t3m
from sage.all import matrix

def relabled_tri():
    M = snappy.Manifold('K13n586')
    tet_perm = [1, 3, 6, 0, 4, 8, 9, 2, 5, 7]
    vertex_perms = [[3, 2, 1, 0],
                    [2, 3, 0, 1],
                    [0, 2, 3, 1],
                    [2, 3, 0, 1],
                    [2, 1, 3, 0],
                    [0, 3, 1, 2],
                    [3, 2, 1, 0],
                    [3, 1, 0, 2],
                    [0, 3, 1, 2],
                    [1, 3, 2, 0]]        
    vertex_perms = [t3m.Perm4(perm) for perm in vertex_perms]
    assert all(p.sign() == 0 for p in vertex_perms)    
    data = M._get_tetrahedra_gluing_data()
    new_data = len(data)*[None]
    for t0, (tets, gluings) in enumerate(data):
        new_tets = 4*[None]
        new_gluings = 4*[None]

        p0 = vertex_perms[t0]
        p0inv = t3m.inv(p0)
        for i, t1 in enumerate(tets):
            new_tets[p0[i]] = tet_perm[t1]

            p1 = vertex_perms[t1]
            curr_perm = t3m.Perm4(gluings[i])
            new_gluings[p0[i]] = (p1 * curr_perm * p0inv).tuple()

        new_data[tet_perm[t0]] = [new_tets, new_gluings]

    T = t3m.Mcomplex(t3m.mcomplex.tets_from_data(new_data))
    assert T.snappy_triangulation() == M.without_hyperbolic_structure()
    return T

def test_tri():
    #T = relabled_tri()
    #S = T.snappy_triangulation()
    S = snappy.Manifold('K13n586_nice.tri')
    master_isosig = snappy.Manifold('K13n586').triangulation_isosig()
    isometry_sig = snappy.Manifold('K13n586').isometry_signature()
    #assert master_isosig == S.triangulation_isosig()
    #assert isometry_sig == master_isosig.split('_')[0]
    
    R = regina.Triangulation3(S._to_string())
    a, b = R.findAllIsomorphisms(R)
    print(b.detail())
    CS = nscomplex.ConnectedSurfaces(S, euler_bound=-6)
    LW = CS.essential_faces_of_normal_polytope()
    G, F = LW.maximal[0].vertex_surfaces
    print()
    print(matrix(10, 7, F.tri_quad_vector)[:,:5])
    print()
    print(matrix(10, 7, G.tri_quad_vector)[:,:5])
    print()
    print(matrix(10, 7, (F + G).tri_quad_vector)[:,:5])


def print_edge_orientations():
    M = snappy.Manifold('K13n586_nice.tri')
    T = t3m.Mcomplex(M)
    for tet in T.Tetrahedra:
        print(tet)
        for i in range(4):
            for j in range(i+1, 4):
                v_i = t3m.ZeroSubsimplices[i]
                v_j = t3m.ZeroSubsimplices[j]
                edge = tet.Class[v_i|v_j]
                orient = edge.orientation_with_respect_to(tet, v_i, v_j)
                line = 4*" " + "E%d%d: %s %d" % (i, j, edge, orient)
                line = line.replace(' (int)', '')
                print(line)
        print()
                
    return T


if __name__ == '__main__':
    test_tri()
