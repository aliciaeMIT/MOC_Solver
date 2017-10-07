from initializetracks import InitializeTracks
from geometry import Geometry

num_azim = 8 #number of azimuthal angles desired
t = 0.5#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3


tracks = InitializeTracks(num_azim, t, w, h, n_p, r)
tracks.getTracks()

#distance = Geometry(w, h, r)

"""
for i in range(tracks.num_azim2):
    for j in range(int(tracks.ntot[i])):
        distance.findIntersection(tracks.startpoint[i][j], tracks.endpoint[i][j])
"""

#tracks.findIntersection()