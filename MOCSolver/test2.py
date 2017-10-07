from initializetracks import InitializeTracks
from geometry import Geometry

num_azim = 8 #number of azimuthal angles desired
t = 0.1#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3
num_rings = 2

segments = Geometry(w, h, r)
segments.createRings(num_rings)

tracks = InitializeTracks(num_azim, t, w, h, n_p, r, num_rings, segments.ring_radius)
tracks.getTracks()
tracks.getStart()
tracks.getEnd()
tracks.getAngularQuadrature()
tracks.getPolarWeight()
tracks.findIntersection()
tracks.plotTracks()



