#Alicia M. Elliott, 22.212 Fall 2017
#Method of Characteristics solver
#2D pincell, fixed isotropic source in fuel
#Square geometry (for comparison with SN)

from initializetracks import InitializeTracks
from flux import FlatSourceRegion, MethodOfCharacteristics, ConvergenceTest
from math import pi


###################################
########## PROBLEM SETUP ##########
###################################

num_azim = 16                  #number of azimuthal angles desired
t = 0.05                        #track spacing desired, cm
h = 1.6                        #height of pincell
w = 1.6                        #width of pincell
r = 0.8/2                        #fuel pin effective radius (half width of square fuel pin)
n_p = 3                         #number of polar divisions; can be 2 or 3
num_iter_max = 100              #maximum number of iterations on flux
tol = 1e-7                     #tolerance for convergence (using L2 Engineering Norm)
fuelgeom = 'square'

#########################################
########## MATERIAL PROPERTIES ##########
#########################################

q_fuel = 10/ 4 * pi                   #constant isotropic source in fuel
q_mod = 0                       #no source in moderator
ndensity_fuel = 2.2e22          #atoms/cc (UO2)
ndensity_mod = 1.0e21           #at/cc (H2O)
sigma_a = 1e-24                 #moderator absorption cross section (cm^2)
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

sigma_t_fuel = 1e5 #sigma_r * ndensity_fuel
sigma_t_mod =  1    #sigma_a * ndensity_mod


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

setup = InitializeTracks(num_azim, t, w, h, n_p, r, fsr, fuelgeom)
setup.getTrackParams()
setup.makeTracks()
setup.getAngularQuadrature()
setup.getPolarWeight()
if fuelgeom == 'square':
    setup.findIntersection()
else:
    setup.findCircleIntersection()
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

