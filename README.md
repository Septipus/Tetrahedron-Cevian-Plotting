# Tetrahedron-Cevian-Plotting
Exploring plotting 4-dimesional data in a 3-dimensional volume (tetrahedron) using cevian averaging

Phase plots display the change in the states of a dynamic system with
respect to themselves.

for example the state of a system with states:
    [s0, s1, s2, ..... sx] with n = timestep number
could be ploted in 2 dimensions by iteratively vaying the state values
along the x & y axes in the following manner:
               x ,    y
    p0      = [s0,    s1]
    p1      = [s1,    s2]
    p2      = [s2,    s3]
            .
            .
            .
    p(n-1)  = [s(n-1), s(n)]
    
The same can be applied to the same data but in 2 dimensions by increasing
the axes count to 3 (x, y, z) and applying the same technique above:
               x ,    y
    p0      = [s0,    s1,     s2]
    p1      = [s1,    s2,     s3]
    p2      = [s2,    s3,     s4]
            .
            .
            .
    p(n-2)  = [s(n-2), s(n-1), s(n)]
    note how the number of points decreases given the # of axes used
    
This however limits the number of axes available for plotting given the
maximum of 3 dimensions available for visualisation with standard
cartesian coordinates.

A re-interpretation of this process in triangular space allows for plotting
of higher dimensional data in the 3 cartesian dimensions available.

In the case of the 3D phase ploting, the data is in the form of:
    p(n-2)  = [s(n-2), s(n-1), s(n)]
    
If the 3 data points are normalised by finding the weight of every state
    w(n-2) = [s(n-2), s(n-1), s(n)] / sum(p(n-2))
    it is evident that 
    sum(w(n-2)) = 1

As such, the values of w(n-2) can be interpreted as the proportional
weights of p(n-2) that can represent the cevian proportions of a point
plotted inside the equilateral triangle given that the weighted cevian 
coordinates of the verticies of a triangle correspond to:
    (1,0,0), (0,1,0), (0,0,1) with the centroid at (1/3, 1/3, 1/3)

Since these three points are representable on the standard (x,y) plane,
for example by considering them as the points of an equilateral triangle
laying on a unit circle, the reinterpreted data is effectively mapped from a
3D space onto a 2D one.

The same holds for the next higher order equivalent. In this case a 4D set
of phase points, can be represented as weighted points within a tetrahedron
and as such are visualisable in standard 3D cartesian coordiantes. In this
case a the circumcribed sphere of a tetrahedron is used with cevian
verticies:
    (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)
    with the centroid at (1/4, 1/4, 1/4, 1/4)

This class transfroms an arbitrary set of generated state values for an 
arbitrary number of states, and plots each of these value sets onto either
a 2D representation of 3 phase state length (in a triangle) or a 3D
representation of 4 phase state length (in a tetrahedron)
