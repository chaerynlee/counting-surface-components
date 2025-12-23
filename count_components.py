import pickle
import snappy
import regina
import pandas as pd
import os
from nscomplex_updated import *
from Orbits import orbits
import itertools

class SurfacetoOrbit:
    '''
    Finds the intersection of a normal surface with the 1, 2 skeleton of the triangulation in which it lies in
    and changes them into an integer interval and collection of pairings
    '''
    def __init__(self, surface):
        self.surface = surface  # this is a regina normal surface
        self.triangulation = surface.triangulation() # this is a regina triangulation
        self.pairings = []
        self.interval = None
        self.interval_divided = []
        self._find_intersections1d()
        self._find_intersections2d()

    def _find_intersections1d(self):
        # list of edge weights of given surface
        edge_weights = [self.surface.edgeWeight(n) for n in range(self.triangulation.countEdges())]
        edge_weights_int = [int(edge_weight.stringValue()) for edge_weight in edge_weights]
        # total number of intersections with edges(1-skeleton)
        num_int = sum(edge_weights_int)
        # single integer interval corresponding to all intersections with edges, starts at 1 not 0
        self.interval = orbits.Interval(1, num_int)

        # list of integer lists where each list contains integers corresponding to the given edge
        # e.g. if edge0 has 3 intersections then we have [[1, 2, 3], [...], ...]
        self.interval_divided = []
        for i, n in enumerate(edge_weights_int):
            if n == 0:
                self.interval_divided.append([])
            else:
                start = sum(edge_weights_int[:i]) + 1
                edge = list(range(start, start + edge_weights_int[i]))
                self.interval_divided.append(edge)

    def oriented_edges(self, triangulation):
        triangles = triangulation.triangles()
        oriented_edges = []
        for tri in triangles:
            tri_info = []  # list of edge and orientation information of given triangle (always has 3 entries)
            for i in range(3):
                edge = tri.edge(i)  # edge index in the triangulation
                mapping = tri.edgeMapping(i)
                # the first two entries correspond to the vertex indices
                # 0: 12 or 21, 1: 02 or 20, 2: 01 or 10 depending on the orientation
                if (mapping[0], mapping[1]) in [(0, 1), (1, 2), (2, 0)]:
                    orientation = 1
                else:
                    orientation = -1
                tri_info.append((int(edge.index()), orientation))
            oriented_edges.append(tri_info)  # list of lists of information for each triangle (has (# of faces) entries which are lists)
        return oriented_edges

    def _find_intersections2d(self):
        ori_edges = self.oriented_edges(self.triangulation)
        triangles = self.triangulation.triangles()

        for t, tri in enumerate(triangles):
            for i in range(3):
                width_regina = self.surface.arcs(tri.index(), i)
                width = int(width_regina.stringValue())
                if width == 0:
                    continue
                else:
                    domain_edge = tri.edge((i + 1) % 3).index()
                    range_edge = tri.edge((i + 2) % 3).index()
                    domain_ori = ori_edges[t][(i + 1) % 3][1]
                    range_ori = ori_edges[t][(i + 2) % 3][1]
                    if domain_ori == 1:
                        domain_interval = [self.interval_divided[domain_edge][-width],
                                           self.interval_divided[domain_edge][-1]]
                    elif domain_ori == -1:
                        domain_interval = [self.interval_divided[domain_edge][0],
                                           self.interval_divided[domain_edge][width - 1]]
                    if range_ori == 1:
                        range_interval = [self.interval_divided[range_edge][0],
                                          self.interval_divided[range_edge][width - 1]]
                    elif range_ori == -1:
                        range_interval = [self.interval_divided[range_edge][-width],
                                          self.interval_divided[range_edge][-1]]
                    if domain_ori == range_ori:
                        self.pairings.append(orbits.Flip(domain_interval, range_interval))
                    else:
                        self.pairings.append(orbits.Shift(domain_interval, range_interval))

    def countcomponents(self):
        if self.pairings == []:
            return self.interval.width
        else:
            G = orbits.Pseudogroup(self.pairings, self.interval)
            return G.reduce()


def main():
    M = snappy.Manifold('t12647')
    CS = connected_surfaces.ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    vertex_surfaces = [face.vertex_surfaces for face in LW.maximal]
    F = vertex_surfaces[0][0]
    G = vertex_surfaces[0][1]
    print(type(F))
    # S = 5 * F + 2 * G
    # SO = SurfacetoOrbit(S.surface)
    # print('interval', SO.interval)
    # print('pairings', len(SO.pairings))
    # for p in sorted(SO.pairings):
    #     print(p)
    # print()
    # print('number of components', SO.countcomponents())

    for comb in itertools.product(range(1, 10), repeat=2):
        S = comb[0] * F + comb[1] * G
        SO = SurfacetoOrbit(S.surface)
        print(comb, 'number of components', SO.countcomponents())

def check(manifold, genus):
    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == manifold].values[0]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])
    T = regina.Triangulation3(TS)
    LWC_info = df.iloc[i, df.columns.get_loc('max_faces')]
    vector_info = df.iloc[i, df.columns.get_loc('vertex_surfaces')]
    genera_info = df.iloc[i, df.columns.get_loc('vertex_genera')]

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

        AF = faces.AdmissibleFace(dim, vertex_surfaces_ns)
        solutions = AF.surfaces_of_potential_genus_in_interior(genus)

        # for sol in solutions:
        #     SO = SurfacetoOrbit(sol)
        #     print('surface', sol)
        #     print(SO.surface.isConnected())
        #     print('number of components', SO.countcomponents())
        #     print()
        S = solutions[3]
        SO = SurfacetoOrbit(S)
        print('surface', S)
        print(SO.surface.isConnected())
        print('number of components', SO.countcomponents())

def check_regina_vs_snappy(manifold):
    df = pd.read_csv(os.getcwd() + '/very_large_combined.csv')
    i = df.index[df['name'] == manifold].values[0]
    TS = snappy.Manifold(df.iloc[i, df.columns.get_loc('tri_used')])

    print('from snappy')
    CS = connected_surfaces.ConnectedSurfaces(TS, -6)
    LW = CS.essential_faces_of_normal_polytope()
    vertex_surfaces = [face.vertex_surfaces for face in LW.maximal]
    print('LW-complex', LW.maximal)
    print('vertex surfaces')
    for F in vertex_surfaces:
        for S in F:
            print(S, S.quad_vector)

    print('from regina')
    T = regina.Triangulation3(TS)
    CS = connected_surfaces.ConnectedSurfaces(T, -6)
    LW = CS.essential_faces_of_normal_polytope()
    vertex_surfaces = [face.vertex_surfaces for face in LW.maximal]
    print('LW-complex', LW.maximal)
    print('vertex surfaces')
    for F in vertex_surfaces:
        for S in F:
            print(S, S.quad_vector)

def test_group():
    I = orbits.Interval(1, 2)
    pairings = [orbits.Pairing(orbits.Interval(1, 1), orbits.Isometry(3, 1))]
    # pairings = []
    G = orbits.Pseudogroup(pairings, I)
    print(G)
    print(G.reduce())


if __name__ == '__main__':
    main()

    # example that didn't work
    # I = orbits.Interval(1, 6)
    # pairings = orbits.Flip((3, 4), (5, 6)), orbits.Shift((2, 4), (4, 6)), orbits.Shift((1, 3), (3, 5))
    # G = orbits.Pseudogroup(pairings, I)
    # print('number of components', G.reduce())

    # p = orbits.Shift((2, 4), (4, 6))
    # print('p', p)
    # q = orbits.Shift((1, 3), (3, 5))
    # print('q', q)
    # g = p.merge(q)
    # print('merged', g)