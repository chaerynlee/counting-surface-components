import snappy
import regina
from nscomplex import *
from Orbits import orbits


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
        #list of edge weights of given surface
        edge_weights = [self.surface.edgeWeight(n) for n in range(self.triangulation.countEdges())]
        edge_weights_int = [int(edge_weight.stringValue()) for edge_weight in edge_weights]
        #total number of intersections with edges(1-skeleton)
        num_int = sum(edge_weights_int)
        #single integer interval corresponding to all intersections with edges, starts at 1 not 0
        self.interval = orbits.Interval(1, num_int)

        #list of integer lists where each list contains integers corresponding to the given edge
        #e.g. if edge0 has 3 intersections then we have [[1, 2, 3], [...], ...]
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
            for i in range(3):
                edge = tri.edge(i)
                mapping = tri.edgeMapping(i)
                if (mapping[0], mapping[1]) in [(0, 1), (1, 2), (2, 0)]:
                    orientation = 1
                else:
                    orientation = -1
                oriented_edges.append((int(edge.index()), orientation))
        return oriented_edges

    def _find_intersections2d(self):
        ori_edges = self.oriented_edges(self.triangulation)
        triangles = self.triangulation.triangles()
        for t, tri in enumerate(triangles):
            for i in range(3):
                width_regina = self.surface.arcs(tri.index(), i)
                width = int(width_regina.stringValue())
                if width == 0:
                    break
                else:
                    domain_edge = tri.edge((i + 1) % 3).index()
                    range_edge = tri.edge((i + 2) % 3).index()
                    domain_ori = ori_edges[domain_edge][1]
                    range_ori = ori_edges[range_edge][1]
                    if domain_ori == 1:
                        domain_interval = [self.interval_divided[domain_edge][-width], self.interval_divided[domain_edge][-1]]
                    elif domain_ori == -1:
                        domain_interval = [self.interval_divided[domain_edge][0], self.interval_divided[domain_edge][width - 1]]
                    if range_ori == 1:
                        range_interval = [self.interval_divided[range_edge][0], self.interval_divided[range_edge][width - 1]]
                    elif range_ori == -1:
                        range_interval = [self.interval_divided[range_edge][-width], self.interval_divided[range_edge][-1]]
                    if domain_ori == range_ori:
                        self.pairings.append(orbits.Flip(domain_interval, range_interval))
                    else:
                        self.pairings.append(orbits.Shift(domain_interval, range_interval))

    def countcomponents(self):
        G = orbits.Pseudogroup(self.pairings, self.interval)
        return G.reduce()


def main():
    import snappy, regina
    import nscomplex
    M = snappy.Manifold('K13n586')
    CS = nscomplex.ConnectedSurfaces(M, -6)
    LW = CS.essential_faces_of_normal_polytope()
    vertex_surfaces = [face.vertex_surfaces for face in LW.maximal]
    F = vertex_surfaces[0][0]
    G = vertex_surfaces[0][1]
    S = 2 * F + 4 * G
    SO = SurfacetoOrbit(S.surface)
    print(SO.pairings, SO.interval)
    print(SO.countcomponents())

if __name__ == '__main__':
    main()