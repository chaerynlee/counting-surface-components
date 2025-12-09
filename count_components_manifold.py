import snappy
import regina
from count_components import *
from nscomplex_updated import *
import orbits_manifold
import pickle
import networkx as nx
import matplotlib.pyplot as plt

class SurfacetoOrbit:
    '''
    Constructs a general linear combination given a list of surfaces (will be vertex surfaces of a certain face in a LW complex)
    Finds the intersection of this linear combination normal surface with the 1, 2 skeleton of the triangulation in which it lies in
    and changes them into an integer interval and collection of pairings
    '''
    def __init__(self, surface_list):
        self.num_vertex = len(surface_list)
        self.R = []  # list of variables as Polynomial objects, named x0, x1, ...
        for n in range(self.num_vertex):
            var_name = 'x' + str(n)
            var_poly = orbits_manifold.Polynomial(0, **{var_name:1})
            self.R.append(var_poly)
        self.vertex_surfaces = surface_list  # this is a list of regina normal surface
        self.triangulation = self.vertex_surfaces[0].triangulation() # this is a regina triangulation
        self.pairings = []
        self.interval = None
        self.interval_divided = []
        # self.edge_connection = nx.Graph()
        # self.edge_connection.add_nodes_from(range(self.triangulation.countEdges()))
        self._find_intersections1d()
        self._find_intersections2d()

    def _find_intersections1d(self):
        # list of edge weights for each vertex surface
        edge_weights_list = []
        for S in self.vertex_surfaces:
            edge_weights_regina = [S.edgeWeight(n) for n in range(self.triangulation.countEdges())]
            edge_weights = [orbits_manifold.Polynomial(int(e.stringValue())) for e in edge_weights_regina]
            edge_weights_list.append(edge_weights)

        edge_weights_int = []
        for n in range(self.triangulation.countEdges()):
            sum_edgeweight = orbits_manifold.Polynomial(0)
            for k in range(self.num_vertex):
                sum_edgeweight += self.R[k] * edge_weights_list[k][n]
            edge_weights_int.append(sum_edgeweight)

        # total number of intersections with edges(1-skeleton)
        total_count = orbits_manifold.Polynomial(0)
        for p in edge_weights_int:
            total_count += p

        # single interval corresponding to all intersections with edges, starts at 0 since it is continuous
        self.interval = orbits_manifold.Interval(orbits_manifold.Polynomial(0), total_count)
        # was self.interval = orbits_manifold.Interval(orbits_manifold.Polynomial(1), total_count)

        # list of integer lists where each list contains the start and end integers corresponding to the given edge
        # e.g. if edge0 has 3 intersections then we have [[0, 3], [...], ...]
        self.interval_divided = []
        for i, n in enumerate(edge_weights_int):
            if n == orbits_manifold.Polynomial(0):
                self.interval_divided.append([])
            else:
                exclude_last = orbits_manifold.Polynomial(0)
                for p in edge_weights_int[:i]:
                    exclude_last += p
                start = exclude_last
                # start = exclude_last + orbits_manifold.Polynomial(1)
                end = start + edge_weights_int[i]
                # end = start + edge_weights_int[i] - orbits_manifold.Polynomial(1)
                self.interval_divided.append([start, end])

    def oriented_edges(self, triangulation):
        triangles = triangulation.triangles()
        oriented_edges = []
        for tri in triangles:
            edge_info = []
            for i in range(3):
                edge = tri.edge(i)
                mapping = tri.edgeMapping(i)
                if (mapping[0], mapping[1]) in [(0, 1), (1, 2), (2, 0)]:
                    orientation = 1
                else:
                    orientation = -1
                edge_info.append((int(edge.index()), orientation))
            oriented_edges.append(edge_info)
        return oriented_edges

    def _find_intersections2d(self):
        ori_edges = self.oriented_edges(self.triangulation)
        triangles = self.triangulation.triangles()

        for t, tri in enumerate(triangles):
            for i in range(3):
                width = orbits_manifold.Polynomial(0)
                for n in range(self.num_vertex):
                    width_regina = self.vertex_surfaces[n].arcs(tri.index(), i)
                    width += self.R[n] * orbits_manifold.Polynomial(int(width_regina.stringValue()))

                if width == orbits_manifold.Polynomial(0):
                    continue
                else:
                    domain_edge = ori_edges[t][(i + 1) % 3][0]
                    range_edge = ori_edges[t][(i + 2) % 3][0]
                    # self.edge_connection.add_edge(domain_edge, range_edge)
                    domain_ori = ori_edges[t][(i + 1) % 3][1]
                    range_ori = ori_edges[t][(i + 2) % 3][1]
                    if domain_ori == 1:
                        domain_interval = orbits_manifold.Interval(self.interval_divided[domain_edge][1] - width, self.interval_divided[domain_edge][1])
                        # domain_interval = orbits_manifold.Interval(self.interval_divided[domain_edge][1] - width + orbits_manifold.Polynomial(1), self.interval_divided[domain_edge][1])
                    elif domain_ori == -1:
                        domain_interval = orbits_manifold.Interval(self.interval_divided[domain_edge][0], self.interval_divided[domain_edge][0] + width)
                        # domain_interval = orbits_manifold.Interval(self.interval_divided[domain_edge][0], self.interval_divided[domain_edge][0] + width - orbits_manifold.Polynomial(1))
                    if range_ori == 1:
                        range_interval = orbits_manifold.Interval(self.interval_divided[range_edge][0], self.interval_divided[range_edge][0] + width)
                        # range_interval = orbits_manifold.Interval(self.interval_divided[range_edge][0], self.interval_divided[range_edge][0] + width - orbits_manifold.Polynomial(1))
                    elif range_ori == -1:
                        range_interval = orbits_manifold.Interval(self.interval_divided[range_edge][1] - width, self.interval_divided[range_edge][1])
                        # range_interval = orbits_manifold.Interval(self.interval_divided[range_edge][1] - width + orbits_manifold.Polynomial(1), self.interval_divided[range_edge][1])
                    if domain_ori == range_ori:
                        self.pairings.append(orbits_manifold.Flip(domain_interval, range_interval, domain_edge, range_edge))
                    else:
                        self.pairings.append(orbits_manifold.Shift(domain_interval, range_interval, domain_edge, range_edge))

    def countcomponents(self):
        G = orbits_manifold.Pseudogroup(self.pairings, self.interval, self.interval_divided)
        return G.reduce()


class SurfaceComponentCount:
    def __init__(self, manifold):
        self.manifold = manifold
        CS = ConnectedSurfaces(self.manifold)
        LW = CS.essential_faces_of_normal_polytope()
        lw_faces = LW.maximal
        self.vertex_sfces = [face.vertex_surfaces for face in lw_faces]

    def surface_in_face(self, v_sfce_list):
        pass


def main():
    import snappy, regina
    import nscomplex_updated
    # M = snappy.Manifold('v2946')
    M = snappy.Manifold('t12647')
    CS = ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    LW_faces = LW.maximal
    test_list = LW_faces[0].vertex_surfaces
    test_regina_list = [S.surface for S in test_list]
    SO = SurfacetoOrbit(test_regina_list)
    print('interval:', SO.interval)

    print('divided interval:')
    for i, subint in enumerate(SO.interval_divided):
        print(i, subint, subint[1] - subint[0])

    # print('edge connection:', SO.edge_connection.edges)
    # nx.draw(SO.edge_connection, with_labels=True)
    # plt.savefig('edgeconnection.png', dpi=300)
    # plt.close()

    print('pairings:')
    for pairing in SO.pairings:
        print(pairing)

    # Other examples: our manifold has the same a_g(M) pattern as the second one in the table in Section8 of DGR
    # look in very_large_combined.csv in nscomplex for manifolds with the same a_g(M) pattern

if __name__ == '__main__':
    main()