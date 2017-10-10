from initializetracks import InitializeTracks
from geometry import Geometry
from flatsource import FlatSourceApproximation
from flux import MOCFlux

num_azim = 4 #number of azimuthal angles desired
t = 0.3#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3
num_rings = 0 #number of rings to create to divide fuel into FSRs
q_fuel = 0 #constant isotropic source (divide by 1/4pi in later code)
q_mod = 0 #no source in moderator


sigma_a = 1
sigma_r = 1 #(11.4 + 8)* 1e24
sigma_t_fuel = sigma_r #fill these in with real numbers later.
sigma_t_mod = sigma_a

psi_in = 0 #initial guess for flux in

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
#tracks.plotFluxPasses()
#tracks.plotTracks()
#tracks.plotSegments()
tracks.getFSRVolumes()
tracks.getFSRAreas()

source = FlatSourceApproximation(q_fuel, q_mod, tracks.segsource)
print "Determining source values for segment regions..."
source.computeSource(tracks.num_segments, sigma_a, sigma_r)

flux = MOCFlux(tracks.num_segments, sigma_t_fuel, sigma_t_mod, tracks.seglength, tracks.tracklengths, tracks.n_p, tracks.num_azim2, tracks.ntot, source.qseg, tracks.segsource)



flux.totalAngularFlux(psi_in, tracks.segsource)
flux.scalarFlux(source.sigma_FSR, tracks.segarea, tracks.omega_m, tracks.omega_p, tracks.t_eff, tracks.sintheta_p, tracks.segangle)
#tracks.plotScalarFlux(flux.psi_scalar)

