from initializetracks import InitializeTracks
from flux import FlatSourceRegion, MethodOfCharacteristics, ConvergenceTest


num_azim = 16 #number of azimuthal angles desired
t = 0.05#track spacing desired, cm
h = 1.26 #height of pincell
w = 1.26 #width of pincell
r = 0.4 #fuel pin radius
n_p = 2 #number of polar divisions; can be 2 or 3
num_rings = 0 #number of rings to create to divide fuel into FSRs
q_fuel = 100 #constant isotropic source
q_mod = 0 #no source in moderator

ndensity_fuel = 2.2e22 #atoms/cc
ndensity_mod = 1.0e21 #at/cc

sigma_a = 1000
sigma_r = (11.4 + 8)* 1e-24

num_iter_max = 1
tol = 1e-10


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

#refined = Geometry(w, h, r, num_rings)
#refined.createRings()
#refined.divideSectors()
#refined.plotFSRs()

###############################################
########## SETUP FLAT SOURCE REGIONS ##########
###############################################

fuel = FlatSourceRegion(q_fuel, sigma_t_fuel)
mod = FlatSourceRegion(q_mod, sigma_t_mod)
fsr = [fuel, mod]

#####################################
########## GENERATE TRACKS ##########
#####################################

setup = InitializeTracks(num_azim, t, w, h, n_p, r, fsr)
setup.getTrackParams()
setup.makeTracks()
setup.getAngularQuadrature()
setup.getPolarWeight()
setup.findIntersection()
setup.reflectRays()
setup.getFSRVolumes(fuel, mod)
#setup.getTrackLinkCoords()


##############################
########## PLOTTING ##########
##############################

#tracks.plotTracks()
#tracks.plotTrackLinking()
#tracks.plotSegments()

######################################
########## SOLVE FOR FLUXES ##########
######################################

flux = MethodOfCharacteristics(sigma_t_fuel, sigma_t_mod, fsr, setup, check)
converged = False
dancoff_flux0 = 0

while not converged:
    flux.solveFlux(num_iter_max)

    psi_in1 = flux.flux_out
    if iter > 1:
        converged = check.isConverged(flux.psi_scalar, psi_old, tol)
    else:
        dancoff_flux0 = flux.flux_out
    psi_old = flux.psi_scalar
dancoff_flux1 = flux.flux_out

if test_dancoff:
    temp = 4 * math.pi * q_fuel / sigma_t_fuel
    tempa = (-1 * temp * 2 / 3) + (5 / 3) * dancoff_flux1
    temp1 = ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux1)
    temp2 = ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux0)
    temp3 = (4 * math.pi - dancoff_flux1) / (4 * math.pi - dancoff_flux0)
    temp4 = temp1 / temp2
    # temp5 = dancoff_flux1/dancoff_flux0
    dancoff_c = 1 - ((4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux1) / (
    (4 * math.pi * q_fuel / sigma_t_fuel) - dancoff_flux0)
    dcheck = 1 - (temp1 / temp2)
    print "Flux 0/flux conv = %f \t/\t %f" % (dancoff_flux0, dancoff_flux1)
    print "Dancoff factor: %f" % (dancoff_c)

print "Iterations to convergence: %d" % (iter)



#tracks.plotScalarFlux(flux.psi_scalar)

