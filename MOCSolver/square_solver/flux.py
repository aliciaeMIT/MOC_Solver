#Alicia M. Elliott, 22.212 Fall 2017
#Method of Characteristics solver
#2D pincell, fixed isotropic source in fuel

import math


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

    def exponentialTerm(self, length, region, p): #calculates exponential term for a given segment and polar angle
        ans = math.exp(-1 * self.segmentXS(region) * length / self.setup.sintheta_p[p])
        return ans

    def preCalculate(self, n_p): #precalculates exponential terms for each segment and polar angle
        for i in range(self.setup.num_azim2):
            for track in self.tracks[i]:
                for s in track.segments:
                    for p in range(n_p):
                        region = s.region
                        if region == 0:  # moderator
                            q_seg = self.regions[1].source
                        elif region == 1:  # fuel
                            q_seg = self.regions[0].source
                        length = s.length
                        s.exponential.append(self.exponentialTerm(length, region, p))


    def angularFlux(self, flux_in, s, p): #s = segment
        region = s.region
        if region == 0: #moderator
            q_seg = self.regions[1].source
        elif region == 1: #fuel
            q_seg = self.regions[0].source
        delta_psi = (flux_in - q_seg / self.segmentXS(region)) * (1 - s.exponential[p])
        return delta_psi

    def segmentXS(self, region):
        if region == 0:
            sigma_t = self.sigma_t_mod
        elif region == 1:
            sigma_t = self.sigma_t_fuel
        return sigma_t

    def solveFlux(self, num_iter_tot, tol):

        stp = self.setup
        self.preCalculate(stp.n_p)
        num_iter = 0
        delta_flux = 0
        fuel = self.regions[0]
        mod = self.regions[1]
        converged = False
        flux_old = []

        print "Solving for fluxes...\n"
        while not converged:
            #initialize scalar flux accumulators
            fuel.flux = 0
            mod.flux = 0

            for p in range(stp.n_p):                                       #loop over polar angles
                for i in range(stp.num_azim2):                             #loop over azimuthal angles
                    for track in self.tracks[i]:                           #loop over all tracks

                        track.quadwt = stp.omega_m[i] * stp.t_eff[i] * stp.omega_p[p] * stp.sintheta_p[p]  * 4 * math.pi
                        flux = track.flux_in[0, p]

                        if num_iter == 0:                                  #initial flux on boundaries is zero
                            flux = 0
                        #print "Forward tracking:"
                        for s in track.segments:                           #loop over segments
                            region = s.region
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux

                            #print "region %g \t flux out %.1e \t\t delta flux %.1e" % (region, flux, delta_flux)
                            if region == 1:
                                fuel.flux += delta_flux * track.quadwt/2
                            elif region == 0:
                                mod.flux += delta_flux * track.quadwt/2
                            else:
                                print "Error in scalar flux calculation (forward)"

                        track.track_out.flux_in[track.refl_out, p] = flux

                        #track in reverse now
                        flux = track.flux_in[1,p]
                        if num_iter == 0:
                            flux = 0
                        #print "Reverse tracking:"
                        for s in track.segments[::-1]:
                            region = s.region
                            delta_flux = self.angularFlux(flux, s, p)
                            flux -= delta_flux

                           # print "region %g \t flux out %.1e \t\t delta flux %.1e" % (region, flux, delta_flux)
                            if region == 1:
                                fuel.flux += delta_flux * track.quadwt/2
                            elif region == 0:
                                mod.flux += delta_flux * track.quadwt/2
                            else:
                                print "Error in scalar flux calculation (reverse)"
                        track.track_in.flux_in[track.refl_in, p] = flux

            fuel.flux = 4 * math.pi * fuel.source / self.segmentXS(1) + fuel.flux / (self.segmentXS(1) * fuel.area)
            mod.flux = 4 * math.pi * mod.source / self.segmentXS(0) + mod.flux / (self.segmentXS(0) * mod.area)

            if num_iter == 0:
                self.dancoff_flux0 = fuel.flux
                print "Dancoff first pass flux: %f \n"%(self.dancoff_flux0)
                num_iter+=1
            elif num_iter >= 1:
                print "Checking convergence for iteration %d" % (num_iter)
                converged = self.check.isConverged([fuel.flux, mod.flux], flux_old, tol)

                if converged:
                    print "\nNumber of iterations to convergence: %d" % (num_iter+1)
                    dancoff = self.check.computeDancoff(self.dancoff_flux0, fuel.flux, fuel.source,  self.segmentXS(1))
                    print "Dancoff factor: %f" %(dancoff)
                elif num_iter >= num_iter_tot:
                    print "Total number of iterations reached before convergence. Break."
                    converged = True
                    dancoff = self.check.computeDancoff(self.dancoff_flux0, fuel.flux, fuel.source, self.segmentXS(1))
                    print "Dancoff factor: %f" % (dancoff)
                else:
                    num_iter +=1

            #clear previous arrays; set new values equal to old values for next run
            flux_old[:] = []
            flux_old = [fuel.flux, mod.flux]

        #normalize fuel flux to 1
        mod.flux /= fuel.flux
        fuel.flux /= fuel.flux

        print "\nSCALAR FLUX\n-----------" \
              "\nFuel = \t\t\t%g \nModerator = \t%g" \
              "\nNumber of iterations: %d" % (fuel.flux, mod.flux, num_iter+1)

        #stp.plotScalarFlux(fuel.flux, mod.flux)


class ConvergenceTest(object):
    def __init__(self):
        """
        class for testing convergence using the L2 engineering norm
        n is the iteration index, i is the vector content index, I is total number of entries in vector.
        """

    def isConverged(self, vec_n, vec_n1, epsilon):
        sum1 = 0
        for i in range(len(vec_n)):
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
        print "This test case should yield a flat flux profile."
        return qfuel, qmod

    def sourceXSConstTest(self, qfuel, sigma_fuel):
        # angular flux everywhere should equal q/sigma
        qmod = qfuel
        sigma_mod = sigma_fuel
        print "Angular flux should equal %f everywhere when converged" %(qmod/sigma_mod)
        return qmod, sigma_mod

    def dancoffFactor(self, qfuel):
        qfuel = qfuel
        qmod = 0
        sigma_fuel = 1e5
        return qfuel, qmod, sigma_fuel

    def computeDancoff(self, phi_1, phi_fin, source, xs):
        const = 4 * math.pi * source / xs
        return 1 - (const - phi_fin)/ (const - phi_1)


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
