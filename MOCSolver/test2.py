from initializetracks import InitializeTracks
from flux import FlatSourceRegion, MethodOfCharacteristics, ConvergenceTest



###################################
########## PROBLEM SETUP ##########
###################################

num_azim = 16                   #number of azimuthal angles desired
t = 0.05                        #track spacing desired, cm
h = 1.26                        #height of pincell
w = 1.26                        #width of pincell
r = 0.4                         #fuel pin radius
n_p = 3                        #number of polar divisions; can be 2 or 3
num_iter_max = 150              #maximum number of iterations on flux
tol = 1e-10                     #tolerance for convergence (using L2 Engineering Norm)

#########################################
########## MATERIAL PROPERTIES ##########
#########################################

q_fuel = 1e5                     #constant isotropic source in fuel
q_mod = 0                       #no source in moderator
ndensity_fuel = 2.2e22          #atoms/cc (UO2)
ndensity_mod = 1.0e21           #at/cc (H2O)
sigma_a = 1e-24#5 * 1e-24                  #moderator absorption cross section (cm^2)
sigma_r = (11.4 + 8)* 1e-24     #fuel absorption cross section (cm^2)


#######################################
########## CHOOSE TEST CASES ##########
#######################################

test_sourcexsconst = False
test_qpropto = False
test_dancoff = False


################################################
########## MACROSCOPIC CROSS SECTIONS ##########
################################################

sigma_t_fuel = 1e5#sigma_r * ndensity_fuel
sigma_t_mod = 1#sigma_a * ndensity_mod


####################################
########## RUN TEST CASES ##########
####################################
check = ConvergenceTest()

if test_sourcexsconst:
    q_mod, sigma_t_mod = check.sourceXSConstTest(q_fuel, sigma_t_fuel)

if test_qpropto:
    q_fuel, q_mod = check.sourceProptoXSTest(sigma_t_fuel, sigma_t_mod)

if test_dancoff:
    q_fuel, q_mod, sigma_t_fuel = check.dancoffFactor(q_fuel)

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

#setup.plotTracks()
#setup.plotTrackLinking()
#setup.plotSegments()

######################################
########## SOLVE FOR FLUXES ##########
######################################

flux = MethodOfCharacteristics(sigma_t_fuel, sigma_t_mod, fsr, setup, check)
flux.solveFlux(num_iter_max, tol)

