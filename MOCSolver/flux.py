import math
import numpy as np

class MOCFlux(object):
    """
    class for calculating angular and scalar flux
    """
    def __init__(self, num_segments, sigma_t_fuel, sigma_t_mod, segment_lengths, track_lengths, n_p, num_azim, ntot, q_seg, seg_source):
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


    def exponentialTerm(self, seg_length, seg_source):
        ans = math.exp(-1 * self.segmentXS(seg_source) * seg_length)
        return ans


    def angularFlux(self, psi_in, q_seg, seg_length, seg_source):

        delta_psi = (psi_in - q_seg / self.segmentXS(seg_source)) * (1 - self.exponentialTerm(seg_length, seg_source))
        return delta_psi


    def segmentXS(self, seg_source):
        sigma_t = 1
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

                        flux_in = flux_in - self.delta_flux[k]
                        print "flux in %f \t delta flux %f \t k %d" %(flux_in, self.delta_flux[k], k)

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
                sum2 = sum1
                sumsum = 1

            self.psi_scalar[k] = ((4 * math.pi) / (sigma[k]))*(self.q_seg[k] + (1 / area[k]) * (sum1))
            scalar1 = self.psi_scalar[k]
            sig1 = sigma[k]
            areaa = area[k]


        max1 = max(self.psi_scalar)
        print "max scalar flux:"
        print max1
        print "\n"

        for k in range(self.num_segments):
            self.psi_scalar[k] = self.psi_scalar[k] / max1
            print " scalar flux %f \t k %d" % (self.psi_scalar[k], k)