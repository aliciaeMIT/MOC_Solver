#Alicia M. Elliott, 22.212 Fall 2017
#Method of Characteristics solver
#2D pincell, fixed isotropic source in fuel
#Square geometry (for comparison with SN)

from initializetracks import InitializeTracks
from flux import MethodOfCharacteristics, ConvergenceTest
from math import pi
import geometry as geom
import plotter
import time


###################################
########## PROBLEM SETUP ##########
###################################

num_azim = 8                  #number of azimuthal angles desired
t = 0.2                        #track spacing desired, cm
h = 1.6                        #height of pincell
w = h                           #width of pincell
r = 0.8/2                        #fuel pin effective radius (half width of square fuel pin)
n_p = 3                         #number of polar divisions; can be 2 or 3
num_iter_max = 100              #maximum number of iterations on flux
tol = 1e-7                     #tolerance for convergence (using L2 Engineering Norm)
fuelgeom = 'square'
pitch = h

spacing = 0.05 #mesh spacing
fwidth = 0.8 #fuel width/height

#########################################
########## MATERIAL PROPERTIES ##########
#########################################

q_fuel = 10/ 4 * pi                   #constant isotropic source in fuel
q_mod = 0                       #no source in moderator
ndensity_fuel = 2.2e22          #atoms/cc (UO2)
ndensity_mod = 1.0e21           #at/cc (H2O)
sigma_a = 1e-24                 #moderator absorption cross section (cm^2)
sigma_r = (11.4 + 8)* 1e-24     #fuel absorption cross section (cm^2)


################################################
########## MACROSCOPIC CROSS SECTIONS ##########
################################################

#total
sigma_fuel_tot = 1e5 #sigma_r * ndensity_fuel
sigma_mod_tot =  1    #sigma_a * ndensity_mod

#absorption
sigma_fuel_abs = 0
sigma_mod_abs = 0

#scatter
sigma_fuel_scatter = sigma_fuel_tot + sigma_fuel_abs
sigma_mod_scatter = sigma_mod_tot + sigma_mod_abs


#########################################
############ RESULTS STORAGE ############
#########################################
#create directory to store plots, results in
timestr = time.strftime("%Y-%m-%d_%H-%M")
pathname = 'plots/' + timestr
plotter.mkdir_p(pathname)
savepath = pathname
resultsfile = pathname + '/' + timestr + '_results'


###############################################
########## SETUP FLAT SOURCE REGIONS ##########
###############################################
#set material objects
fuelmat = geom.Material('fuel', q_fuel, sigma_fuel_tot, sigma_fuel_scatter, sigma_fuel_abs)
moderator = geom.Material('moderator', q_mod, sigma_mod_tot, sigma_mod_scatter, sigma_mod_abs)


fuel = geom.FlatSourceRegion(q_fuel, sigma_fuel_tot)
mod = geom.FlatSourceRegion(q_mod, sigma_mod_tot)
fsr = [fuel, mod]


# setup mesh cells
mesh = geom.Geometry(pitch, spacing, fwidth, fuelmat, moderator)
mesh.setMesh()

cell_width = mesh.getWidth(pitch, spacing)
fuel_width = mesh.getWidth(fwidth, spacing)
plot_cells = mesh.getPlotCells(cell_width, fuel_width)
plotter.plotMaterial(mesh, spacing, plot_cells, savepath)

#####################################
########## GENERATE TRACKS ##########
#####################################
check = ConvergenceTest()
setup = InitializeTracks(num_azim, t, w, h, n_p, r, fsr, fuelgeom)
setup.getTrackParams()
setup.makeTracks()
setup.findAllTrackCellIntersect(mesh.cells, spacing)
setup.plotCellSegments(spacing, savepath)

setup.getAngularQuadrature()
setup.getPolarWeight()
if fuelgeom == 'square':
    setup.findIntersection()
#else:
#    setup.findCircleIntersection()
setup.reflectRays()
setup.getFSRVolumes(fuel, mod, mesh)
#setup.getTrackLinkCoords()

##############################
########## PLOTTING ##########
##############################

#setup.plotTracks(savepath)
#setup.plotTrackLinking(savepath)
#setup.plotSegments()

######################################
########## SOLVE FOR FLUXES ##########
######################################

flux = MethodOfCharacteristics(sigma_fuel_tot, sigma_mod_tot, fsr, setup, check, mesh)
flux.solveFlux(num_iter_max, tol)


midpt = mesh.n_cells/2 - 1
plotter.plotScalarFlux(mesh, 0, mesh.mesh, 0, savepath)
plotter.plotCenterFlux(mesh, mesh.cells, midpt, 0, 0, savepath)
plotter.plotCenterFluxY(mesh, mesh.cells, midpt, 0, 0, savepath)
