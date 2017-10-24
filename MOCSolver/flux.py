import math
import numpy as np

class MOCFlux(object):
    """
    class for calculating angular and scalar flux
    """
    def __init__(self, ):
        self.num_segments = num_segments
        self.psi_angular = np.zeros(2 * num_segments)
        self.psi_scalar = np.zeros(num_segments)
        self.sigma_t_fuel = sigma_t_fuel
        self.sigma_t_mod = sigma_t_mod
        self.seglength = segment_lengths
        self.tracklength = track_lengths
        self.n_p = n_p
        self.num_azim = num_azim
        self.ntot = ntot
        self.q_seg = q_seg
        self.segsource = seg_source
        self.delta_flux = np.zeros((num_segments))
        self.flux_out = 0


    def exponentialTerm(self, seg_length, seg_source):
        ans = math.exp(-1 * self.segmentXS(seg_source) * seg_length)
        return ans


    def angularFlux(self, psi_in, q_seg, seg_length, seg_source):

        delta_psi = (psi_in - q_seg / self.segmentXS(seg_source)) * (1 - self.exponentialTerm(seg_length, seg_source))
        #delta_psi = (q_seg / self.segmentXS(seg_source) - psi_in) * (1 - self.exponentialTerm(seg_length, seg_source))
        return delta_psi


    def segmentXS(self, seg_source):

        if seg_source == 0:
            sigma_t = self.sigma_t_mod
        elif seg_source == 1:
            sigma_t = self.sigma_t_fuel
        return sigma_t

    def totalAngularFlux(self, psi_in, seg_source):

        flux_in = psi_in
        l=0

        i0 = 0
        j0 = 0
        k= 0
        diff =1

        for p in range(self.n_p):                                       #loop over polar angles
            for i in range(i0, self.num_azim):
                #loop over azimuthal angles
                for j in range(j0, int(self.ntot[i])):                           #loop over tracks
                    #for l in [0,1]:                                     #loop forward and backward for each track
                    seglen_count = 0
                    counter = self.tracklength[i][j]
                    count = 0
                    diff = round((counter - seglen_count), 4)

                    while not diff == 0:                 #loop over segments in each track
                        count += 1

                        seglen_count += self.seglength[k]

                        diff = round((counter - seglen_count), 4)

                        self.delta_flux[k] = self.angularFlux(flux_in, self.q_seg[k], self.seglength[k], self.segsource[k])


                        print "flux in %f \t delta flux %f \t k %d" %(flux_in, self.delta_flux[k], k)

                        flux_in = flux_in - self.delta_flux[k]
                        if self.segsource[k] == 1:
                            fout = flux_in
                            self.flux_out = flux_in

                        if count > 1 and diff == 0 and l == 0:
                            l += 1
                            seglen_count = 0
                            k -= 1
                            diff = 1

                        if count > 3 and diff !=0:
                            k -= 2
                        if count > 3 and diff == 0 and l>0:
                            k += 3
                            l=0
                            seglen_count = 0
                            continue
                        else:
                            k += 1
                    if k >= self.num_segments:
                        break
                if k >= self.num_segments:
                    break
            if k >= self.num_segments:

                break


    def scalarFlux(self, sigma, area, omega_m, omega_p, omega_k, sinthetap, segangle):
        for k in range(self.num_segments):
            sum1 = 0
            for p in range(self.n_p):
                index = segangle[k][0]
                sum1 += omega_m[index] * omega_p[p] * omega_k[index] * sinthetap[p] * self.delta_flux[k]
                #for debugging:
                sum2 = sum1


            self.psi_scalar[k] = ((4 * math.pi) / (sigma[k]))*(self.q_seg[k] + (1 / area[k]) * (sum1))
            #below is for debugging
            scalar1 = self.psi_scalar[k]
            sig1 = sigma[k]
            areaa = area[k]


        max1 = max(self.psi_scalar)
        print "max scalar flux:"
        print max1
        print "\n"

        if max1 == 'nan' or round(max1,5) == 0.00000:
            pass
        else:
            for k in range(self.num_segments):
                self.psi_scalar[k] = self.psi_scalar[k] / max1
                print " scalar flux %f \t k %d" % (self.psi_scalar[k], k)



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
    def __init__(self, source):

        """
        This class generates the source for each flat source (and updates it, if not constant isotropic)
        """
        self.volume = 0
        self.flux = 0
        self.source = source
        self.segments = []

    def addRegion(self):
        pass

