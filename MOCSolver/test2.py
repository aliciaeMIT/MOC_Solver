from initializetracks import InitializeTracks
from geometry import Geometry
from flatsource import FlatSourceApproximation
from flux import MOCFlux
from convergence import ConvergenceTest

num_azim = 8 #number of azimuthal angles desired
t = 0.1#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3
num_rings = 0 #number of rings to create to divide fuel into FSRs
q_fuel = 1 #constant isotropic source (divide by 1/4pi in later code)
q_mod = 0 #no source in moderator

ndensity_fuel = 2.2e22 #atoms/cc
ndensity_mod = 1.0e21 #at/cc

sigma_a = 1 * 1e-24
sigma_r = (11.4 + 8)* 1e-24

sigma_t_fuel = 1e5
#sigma_t_fuel = sigma_r * ndensity_fuel
sigma_t_mod = sigma_a * ndensity_mod


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
source.computeSource(tracks.num_segments, sigma_t_mod, sigma_t_fuel)

flux = MOCFlux(tracks.num_segments, sigma_t_fuel, sigma_t_mod, tracks.seglength, tracks.tracklengths, tracks.n_p, tracks.num_azim2, tracks.ntot, source.qseg, tracks.segsource)
psi_in1 = psi_in



tol = 1e-5
#for count in range(30):
while not converged:
    flux.totalAngularFlux(psi_in1, tracks.segsource)
    flux.scalarFlux(source.sigma_FSR, tracks.segarea, tracks.omega_m, tracks.omega_p, tracks.t_eff, tracks.sintheta_p, tracks.segangle)
    psi_in1 = flux.flux_out
    converged = ConvergenceTest.isConverged(flux.psi_scalar, psi_old, tol)
    psi_old = flux.psi_scalar

#tracks.plotScalarFlux(flux.psi_scalar)

