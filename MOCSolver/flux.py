import math
import numpy as np

class MethodOfCharacteristics(object):
    """
    class for calculating angular and scalar flux
    """
    def __init__(self, sigma_t_fuel, sigma_t_mod, regions, setup, check):

        self.sigma_t_fuel = sigma_t_fuel
        self.sigma_t_mod = sigma_t_mod
        self.regions = regions
        self.setup = setup
        self.tracks = setup.tracks
        self.check = check
        #get n_p by self.setup.n_p, numazim2 by self.setup.num_azim2
        #get segments in a track by self.tracks.segments
        #get segments in a region by self.regions.segments

        """
        self.psi_angular = np.zeros(2 * num_segments)
        self.psi_scalar = np.zeros(num_segments)

        segment length = self.setup.tracks.segments[].length        
        
        self.tracklength = track_lengths
        self.n_p = n_p
        self.num_azim = num_azim
        self.ntot = ntot
        self.q_seg = q_seg
        self.segsource = seg_source
        self.delta_flux = np.zeros((num_segments))
        self.flux_out = 0
        """

    def exponentialTerm(self, length, region, p):
        ans = math.exp(-1 * self.segmentXS(region) * length / self.setup.sintheta_p[p])
        return ans


    def angularFlux(self, flux_in, s, p): #s = segment
        region = s.region
        if region == 0: #moderator
            q_seg = self.regions[0].source
        elif region == 1: #fuel
            q_seg = self.regions[1].source
        #q_seg = region.source
        length = s.length
        delta_psi = (flux_in - q_seg / self.segmentXS(region)) * (1 - self.exponentialTerm(length, region, p))
        return delta_psi


    def segmentXS(self, region):
        if region == 0:
            sigma_t = self.sigma_t_mod
        elif region == 1:
            sigma_t = self.sigma_t_fuel
        return sigma_t

    def solveFlux(self, num_iter_tot):
        stp = self.setup
        num_iter = 0
        delta_flux = 0
        fuel = self.regions[1]
        mod = self.regions[0]
        converged = False

        while not converged:                              #num_iter < num_iter_tot:
            for p in range(stp.n_p):                                       #loop over polar angles
                for i in range(stp.num_azim2):                             #loop over azimuthal angles
                    for track in self.tracks[i]:                           #loop over all tracks

                        track.quadwt = stp.omega_m[i] * stp.t_eff[i] * stp.omega_p[p] * stp.sintheta_p[p]
                        flux = track.flux_in[0, p]

                        if num_iter == 0:                                  #initial flux on boundaries is zero
                            flux = 0

                        for s in track.segments:                           #loop over segments
                            region = s.region
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux
                            print "region %g \t flux in %.1e \t\t delta flux %.1e" % (region, flux, delta_flux)
                            if region == 1:
                                fuel.flux += delta_flux * track.quadwt
                            elif region == 0:
                                mod.flux += delta_flux * track.quadwt
                            else:
                                print "Error in scalar flux calculation (forward)"

                        track.track_out.flux_in[track.refl_out, p] = flux

                        #track in reverse now
                        flux = track.flux_in[1,p]
                        if num_iter == 0:
                            flux = 0

                        for s in track.segments[::-1]:
                            region = s.region
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux
                            if region == 1:
                                fuel.flux += delta_flux * track.quadwt
                            elif region == 0:
                                mod.flux += delta_flux * track.quadwt
                            else:
                                print "Error in scalar flux calculation (reverse)"
                        track.track_in.flux_in[track.refl_in, p] = flux

                       # print "flux in %g \t delta flux %g" % (flux, delta_flux)
            if num_iter == 0:
                self.dancoff_flux0 = fuel.flux
            num_iter += 1

        fuel.flux /= fuel.area
        fuel.flux += fuel.source
        fuel.flux *= (4*math.pi)/self.segmentXS(1)


        mod.flux /= mod.area
        mod.flux += mod.source
        mod.flux *= (4*math.pi)/self.segmentXS(0)

        print "\nSCALAR FLUX\n-----------" \
              "\nFuel = \t\t\t%g \nModerator = \t%g" %(fuel.flux, mod.flux)


class ConvergenceTest(object):
    def __init__(self):
        """
        class for testing convergence using the L2 engineering norm
        n is the iteration index, i is the vector content index, I is total number of entries in vector.
        """


    def isConverged(self, vec_n, vec_n1, epsilon):
        sum1 = 0
        print "Convergence error:"
        for i in range(len(vec_n)):
            if round(vec_n[i],7) == 0:
                vec_n[i] = 0.000001
            error1 = ((vec_n[i] - vec_n1[i])/vec_n[i]) ** 2
            print error1
            sum1 += error1
        I = len(vec_n)
        l2 = math.sqrt((1 / I) * sum1)

        if l2 < epsilon:
            print "Converged! l2 %f" %(l2)
            return True
        else:
            print "Not converged; l2 %f" %(l2)
            return False


    def sourceProptoXSTest(self, xsfuel, xsmod):
        """test that sets source in each region proportional to the cross section in that region
        should yield a flat flux profile if code is functioning correctly"""
        qfuel = 2 * xsfuel
        qmod = 2 * xsmod
        return qfuel, qmod

    def sourceXSConstTest(self, qfuel, sigma_fuel):
        # angular flux everywhere should equal q/sigma
        qmod = qfuel
        sigma_mod = sigma_fuel
        print "Angular flux should equal %f everywhere when converged" %(qmod/sigma_mod)
        return qmod, sigma_mod

    def dancoffFactor(self):
        qfuel = 1e3
        qmod = 0
        psi_in = 0
        sigma_fuel = 1e5
        return qfuel, qmod, psi_in, sigma_fuel


class FlatSourceRegion(object):
    def __init__(self, source, xs):

        """
        This class generates the source for each flat source (and updates it, if not constant isotropic)
        """
        self.volume = 0
        self.area = 0
        self.flux = 0
        self.source = source
        self.sigma = xs
        #self.segments = []

    def dancoffFlux(self, dflux):
        pass
