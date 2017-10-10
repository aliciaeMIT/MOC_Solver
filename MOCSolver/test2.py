from initializetracks import InitializeTracks
from geometry import Geometry
from flatsource import FlatSourceApproximation
from flux import MOCFlux

num_azim = 4 #number of azimuthal angles desired
t = 0.1#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 3 #number of polar divisions; can be 2 or 3
num_rings = 0 #number of rings to create to divide fuel into FSRs
q_fuel = 1 #constant isotropic source (divide by 1/4pi in later code)
q_mod = 0 #no source in moderator

sigma_t_fuel = 1 #fill these in with real numbers later.
sigma_t_mod = 1

refined = Geometry(w, h, r, num_rings)
#segments.createRings()
#segments.divideSectors()
#segments.plotFSRs()

tracks = InitializeTracks(num_azim, t, w, h, n_p, r, num_rings, refined.ring_radius)

tracks.getTracks()
tracks.getStart()
tracks.getEnd()
tracks.getAngularQuadrature()
tracks.getPolarWeight()
tracks.findIntersection()
#tracks.plotTracks()
#tracks.plotSegments()
tracks.getFSRVolumes()

source = FlatSourceApproximation(q_fuel, q_mod, tracks.segsource)
print "Determining source values for segment regions..."
source.computeSource(tracks.num_segments)

flux = MOCFlux(tracks.num_segments, sigma_t_fuel, sigma_t_mod, tracks.seglength, tracks.tracklengths, tracks.n_p, tracks.num_azim2, tracks.ntot, source.qseg, tracks.segsource)
k=1
psi_in = 1 #initial guess for flux in
flux.totalAngularFlux(psi_in, tracks.seglength)
"""
for p in range(n_p):
    for i in range(tracks.num_azim2):
        for k in range(tracks.num_segments):

            flux.angularFlux(psi_in, source.qseg[k], tracks.seglength[k], tracks.segsource[k])
"""