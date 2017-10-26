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
        self.exponential = []
        #get n_p by self.setup.n_p, numazim2 by self.setup.num_azim2
        #get segments in a track by self.tracks.segments
        #get segments in a region by self.regions.segments

    def exponentialTerm(self, length, region, p):
        ans = math.exp(-1 * self.segmentXS(region) * length )#/ self.setup.sintheta_p[p])
        return ans

    def preCalculate(self):
        p=1
        for i in range(self.setup.num_azim2):
            for track in self.tracks[i]:
                for s in track.segments:
                    region = s.region
                    if region == 0:  # moderator
                        q_seg = self.regions[0].source
                    elif region == 1:  # fuel
                        q_seg = self.regions[1].source
                    length = s.length
                    exp = self.exponentialTerm(length, region, p)
                    s.exponential = exp
                    #self.exponential.append(0)

    def angularFlux(self, flux_in, s, p): #s = segment
        region = s.region
        if region == 0: #moderator
            q_seg = self.regions[0].source
        elif region == 1: #fuel
            q_seg = self.regions[1].source
        #q_seg = region.source
        length = s.length
        #delta_psi = (flux_in - q_seg / self.segmentXS(region)) * (1 - self.exponentialTerm(length, region, p))
        delta_psi = (flux_in - q_seg / self.segmentXS(region)) * (1 - s.exponential)
        return delta_psi

    def segmentXS(self, region):
        if region == 0:
            sigma_t = self.sigma_t_mod
        elif region == 1:
            sigma_t = self.sigma_t_fuel
        return sigma_t

    def solveFlux(self, num_iter_tot, tol):
        self.preCalculate()
        stp = self.setup
        num_iter = 0
        delta_flux = 0
        fuel = self.regions[1]
        mod = self.regions[0]
        converged = False
        flux_old = []
        flux_tot_old = []
        flux_curr = []
        flux_tot_curr = []


        while not converged:   #num_iter < num_iter_tot:
            for p in range(stp.n_p):                                       #loop over polar angles
                for i in range(stp.num_azim2):                             #loop over azimuthal angles
                    for track in self.tracks[i]:                           #loop over all tracks
                        fuel_old = fuel.flux
                        mod_old = mod.flux

                        #del flux_curr[:]

                        track.quadwt = 4*math.pi* stp.omega_m[i] * stp.t_eff[i] * stp.omega_p[p] * stp.sintheta_p[p]
                        flux = track.flux_in[0, p]

                        if num_iter == 0:                                  #initial flux on boundaries is zero
                            flux = 0

                        for s in track.segments:                           #loop over segments

                            region = s.region
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux

                            #print "region %g \t flux in %.1e \t\t delta flux %.1e" % (region, flux, delta_flux)
                            if region == 1:
                                fuel.flux += delta_flux * track.quadwt
                                flux_curr.append(fuel.flux)
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
                                flux_curr.append(fuel.flux)
                            elif region == 0:
                                mod.flux += delta_flux * track.quadwt
                            else:
                                print "Error in scalar flux calculation (reverse)"

                        track.track_in.flux_in[track.refl_in, p] = flux

            fuel.flux /= fuel.area
            fuel.flux += fuel.source
            fuel.flux *= (4 * math.pi) / self.segmentXS(1)

            mod.flux /= mod.area
            mod.flux += mod.source
            mod.flux *= (4 * math.pi) / self.segmentXS(0)

            if num_iter == 0:
                self.dancoff_flux0 = fuel.flux
                num_iter+=1
            elif num_iter >= 1:
                print "Checking convergence for iteration %d" % (num_iter+1)
                converged = self.check.isConverged(flux_curr, flux_old, tol)

                if converged:
                    print "\nNumber of iterations to convergence: %d" % (num_iter+1)
                    print self.check.computeDancoff(self.dancoff_flux0, fuel.flux)
                elif num_iter >= num_iter_tot:
                    print "Total number of iterations reached before convergence. Break."
                    converged = True
                else:
                    num_iter +=1

            flux_old[:] = []
            flux_old = flux_curr[:]
            flux_curr[:] = []

        print "\nSCALAR FLUX\n-----------" \
              "\nFuel = \t\t\t%g \nModerator = \t%g" \
              "\nNumber of iterations: %d" % (fuel.flux, mod.flux, num_iter+1)


class ConvergenceTest(object):
    def __init__(self):
        """
        class for testing convergence using the L2 engineering norm
        n is the iteration index, i is the vector content index, I is total number of entries in vector.
        """

    def isConverged(self, vec_n, vec_n1, epsilon):
        sum1 = 0
        for i in range(len(vec_n)):
           # if round(vec_n[i],10) == 0:
            #    vec_n[i] = 0.000000000000001
            error1 = ((vec_n[i] - vec_n1[i])/vec_n[i]) ** 2
            sum1 += error1
        I = len(vec_n)
        l2 = math.sqrt(sum1 / I)

        if l2 < epsilon:
            print "Converged! l2 %g" %(l2)
            return True
        else:
            print "Not converged; l2 %g" %(l2)
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
        sigma_fuel = 1e5
        return qfuel, qmod, sigma_fuel

    def computeDancoff(self, phi_1, phi_fin):
        return 1 - phi_fin/phi_1


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
