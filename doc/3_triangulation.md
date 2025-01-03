## 3. Delauney triangulation

Perform Delaunay triangulation on the system so that cavity analysis can be performed. Other methods: Gaussian density analysis, brutal force cutoff analysis.

1. Delauney triangulation can be achieved in two distinct ways: 1. from an initial small unit tetrahedron, insert point from outside of the convex hull once a time untill all points are inserted. 2. from a big enough initial tetrahedron that include all points inside it, insert point from inside of the tetrahedron once a time untill all points are inserted.

2. At PBC, atoms at outlayer need to be dealed seperately.

3. For a polygon with l edges, it is formed by (l-2) triangles. Adding a new point and link the edges to the point, we will have l triangles.
