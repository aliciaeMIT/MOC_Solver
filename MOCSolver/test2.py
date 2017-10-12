from initializetracks import InitializeTracks
from geometry import Geometry
from flatsource import FlatSourceApproximation
from flux import MOCFlux
from convergence import ConvergenceTest
import math

num_azim = 16 #number of azimuthal angles desired
t = 0.05#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 3 #number of polar divisions; can be 2 or 3
num_rings = 0 #number of rings to create to divide fuel into FSRs
q_fuel = 10 #constant isotropic source
q_mod = 0 #no source in moderator

ndensity_fuel = 2.2e22 #atoms/cc
ndensity_mod = 1.0e21 #at/cc

sigma_a = 1000
sigma_r = (11.4 + 8)* 1e-24

psi_in = 0 #initial guess for flux in

check = ConvergenceTest()

#turn test cases on and off
test_sourcexsconst = False
test_qpropto = False
test_dancoff = True


sigma_t_fuel = sigma_r * ndensity_fuel
sigma_t_mod = sigma_a * ndensity_mod


if test_sourcexsconst:
    q_mod, sigma_t_mod = check.sourceXSConstTest(q_fuel, sigma_t_fuel)

if test_qpropto:
    q_fuel, q_mod = check.sourceProptoXSTest(sigma_t_fuel, sigma_t_mod)

if test_dancoff:
    q_fuel, q_mod, psi_in, sigma_t_fuel = check.dancoffFactor()
    #q_fuel = 1
    #q_mod = 0
    #psi_in = 0
    #sigma_t_fuel = 1e5



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

source.computeSource(tracks.num_segments, sigma_t_mod, sigma_t_fuel)

flux = MOCFlux(tracks.num_segments, sigma_t_fuel, sigma_t_mod, tracks.seglength, tracks.tracklengths, tracks.n_p, tracks.num_azim2, tracks.ntot, source.qseg, tracks.segsource)
psi_in1 = psi_in



tol = 1e-10
converged = False
psi_old = flux.psi_scalar

iter = 0

#for count in range(30):
while not converged:
    iter += 1
    flux.totalAngularFlux(psi_in1, tracks.segsource)
    flux.scalarFlux(source.sigma_FSR, tracks.segarea, tracks.omega_m, tracks.omega_p, tracks.t_eff, tracks.sintheta_p, tracks.segangle)
    psi_in1 = flux.flux_out
    if iter > 1:
        converged = check.isConverged(flux.psi_scalar, psi_old, tol)
    else:
        dancoff_flux0 = 0#flux.flux_out
    psi_old = flux.psi_scalar
dancoff_flux1 = flux.flux_out

if test_dancoff:
    temp = 4 * math.pi * q_fuel / sigma_t_fuel
    tempa = (-1*temp * 2 / 3) + (5 / 3 ) * dancoff_flux1
    temp1 = ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux1)
    temp2 = ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux0)
    temp3 = (4 * math.pi - dancoff_flux1)/(4 * math.pi -dancoff_flux0)
    temp4 = temp1 / temp2
    #temp5 = dancoff_flux1/dancoff_flux0
    dancoff_c = 1 - ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux1)/((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux0)
    dcheck = 1 - (temp1/temp2)
    print "Dancoff factor: %f" %(dancoff_c)

print "Iterations to convergence: %d" %(iter)
#tracks.plotScalarFlux(flux.psi_scalar)

