# how to start?ㅠㅠ

M = snappy.Manifold('K13n586')
T = snappy.snap.t3mlite.Mcomplex(M)
T.info()

# tet0    [tet1, tet2, tet3, tet4]
#         [(0, 1, 3, 2), (0, 1, 3, 2), (0, 1, 3, 2), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e0 (int)     E02 : e1 (int)     E12 : e2 (int)
#                E03 : e3 (int)     E31 : e4 (int)     E23 : e4 (int)
# tet1    [tet0, tet5, tet7, tet6]
#         [(0, 1, 3, 2), (0, 1, 3, 2), (0, 1, 3, 2), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e5 (int)     E02 : e6 (int)     E12 : e4 (int)
#                E03 : e1 (int)     E31 : e2 (int)     E23 : e4 (int)
# tet2    [tet6, tet0, tet5, tet8]
#         [(0, 2, 1, 3), (0, 1, 3, 2), (0, 1, 3, 2), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e1 (int)     E02 : e3 (int)     E12 : e7 (int)
#                E03 : e1 (int)     E31 : e8 (int)     E23 : e4 (int)
# tet3    [tet9, tet8, tet7, tet0]
#         [(0, 1, 3, 2), (0, 1, 3, 2), (0, 3, 2, 1), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e0 (int)     E02 : e3 (int)     E12 : e4 (int)
#                E03 : e5 (int)     E31 : e7 (int)     E23 : e9 (int)
# tet4    [tet6, tet5, tet0, tet9]
#         [(1, 2, 3, 0), (1, 2, 3, 0), (0, 1, 3, 2), (2, 0, 3, 1)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e0 (int)     E02 : e9 (int)     E12 : e8 (int)
#                E03 : e1 (int)     E31 : e2 (int)     E23 : e6 (int)
# tet5    [tet9, tet1, tet4, tet2]
#         [(2, 0, 3, 1), (0, 1, 3, 2), (3, 0, 1, 2), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e1 (int)     E02 : e1 (int)     E12 : e8 (int)
#                E03 : e6 (int)     E31 : e9 (int)     E23 : e4 (int)
# tet6    [tet2, tet4, tet1, tet7]
#         [(0, 2, 1, 3), (3, 0, 1, 2), (0, 1, 3, 2), (2, 3, 1, 0)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e5 (int)     E02 : e2 (int)     E12 : e7 (int)
#                E03 : e6 (int)     E31 : e4 (int)     E23 : e8 (int)
# tet7    [tet6, tet8, tet3, tet1]
#         [(3, 2, 0, 1), (1, 3, 0, 2), (0, 3, 2, 1), (0, 1, 3, 2)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e5 (int)     E02 : e1 (int)     E12 : e2 (int)
#                E03 : e0 (int)     E31 : e7 (int)     E23 : e5 (int)
# tet8    [tet9, tet3, tet2, tet7]
#         [(3, 2, 0, 1), (0, 1, 3, 2), (0, 1, 3, 2), (2, 0, 3, 1)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e1 (int)     E02 : e5 (int)     E12 : e0 (int)
#                E03 : e3 (int)     E31 : e7 (int)     E23 : e9 (int)
# tet9    [tet3, tet4, tet5, tet8]
#         [(0, 1, 3, 2), (1, 3, 0, 2), (1, 3, 0, 2), (2, 3, 1, 0)]
#         Vertices: v0 (int) v0 (int) v0 (int) v0 (int)
#         Edges: E01 : e9 (int)     E02 : e0 (int)     E12 : e7 (int)
#                E03 : e8 (int)     E31 : e4 (int)     E23 : e9 (int)
#
# Edges:
# e0 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E01 | F2 | tet0 >  < E01 | F2 | tet3 >  < E03 | F1 | tet7 >
#         < E12 | F0 | tet8 >  < E02 | F1 | tet9 >  < E01 | F2 | tet4 >
# e1 (int)         Edge of valence 9      Endpoints [v0 (int) , v0 (int) ]
#         < E02 | F3 | tet0 >  < E03 | F1 | tet4 >  < E01 | F3 | tet5 >
#         < E01 | F3 | tet2 >  < E01 | F3 | tet8 >  < E02 | F3 | tet7 >
#         < E03 | F1 | tet1 >  < E02 | F3 | tet5 >  < E03 | F1 | tet2 >
# e2 (int)         Edge of valence 5      Endpoints [v0 (int) , v0 (int) ]
#         < E12 | F3 | tet0 >  < E31 | F0 | tet4 >  < E02 | F3 | tet6 >
#         < E12 | F3 | tet7 >  < E31 | F0 | tet1 >
# e3 (int)         Edge of valence 4      Endpoints [v0 (int) , v0 (int) ]
#         < E03 | F1 | tet0 >  < E02 | F3 | tet2 >  < E03 | F1 | tet8 >
#         < E02 | F3 | tet3 >
# e4 (int)         Edge of valence 9      Endpoints [v0 (int) , v0 (int) ]
#         < E31 | F2 | tet0 >  < E12 | F0 | tet3 >  < E31 | F2 | tet9 >
#         < E23 | F1 | tet5 >  < E23 | F0 | tet1 >  < E23 | F1 | tet0 >
#         < E23 | F0 | tet2 >  < E31 | F2 | tet6 >  < E12 | F0 | tet1 >
# e5 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E01 | F2 | tet1 >  < E01 | F2 | tet7 >  < E03 | F1 | tet3 >
#         < E02 | F3 | tet8 >  < E23 | F0 | tet7 >  < E01 | F2 | tet6 >
# e6 (int)         Edge of valence 4      Endpoints [v0 (int) , v0 (int) ]
#         < E02 | F3 | tet1 >  < E03 | F1 | tet6 >  < E23 | F1 | tet4 >
#         < E03 | F1 | tet5 >
# e7 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E12 | F3 | tet2 >  < E31 | F0 | tet8 >  < E12 | F0 | tet9 >
#         < E31 | F2 | tet3 >  < E31 | F0 | tet7 >  < E12 | F0 | tet6 >
# e8 (int)         Edge of valence 5      Endpoints [v0 (int) , v0 (int) ]
#         < E31 | F2 | tet2 >  < E12 | F0 | tet5 >  < E03 | F1 | tet9 >
#         < E12 | F0 | tet4 >  < E23 | F0 | tet6 >
# e9 (int)         Edge of valence 6      Endpoints [v0 (int) , v0 (int) ]
#         < E23 | F1 | tet3 >  < E23 | F0 | tet8 >  < E01 | F2 | tet9 >
#         < E31 | F2 | tet5 >  < E02 | F3 | tet4 >  < E23 | F0 | tet9 >

CS = ConnectedSurfaces(M, -8)