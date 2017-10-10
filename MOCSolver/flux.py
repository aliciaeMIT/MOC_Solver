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
        k = 0
        flux_in = psi_in
        l=0
        #count = 0
        for p in range(self.n_p):                                       #loop over polar angles
            for i in range(self.num_azim):                              #loop over azimuthal angles
                for j in range(int(self.ntot[i])):                           #loop over tracks
                    #for l in [0,1]:                                     #loop forward and backward for each track
                    seglen_count = 0
                    counter = self.tracklength[i][j]
                    count = 0
                    while seglen_count < self.tracklength[i][j]: #loop over segments in each track
                        count += 1

                        seglen_count += self.seglength[k]

                        #delta_flux = self.angularFlux(flux_in, self.q_seg[k], self.seglength[k], self.segsource[k])
                        delta_flux = self.angularFlux(flux_in, 1,1,1)
                        flux_in = flux_in - delta_flux

                        if count > 1 and seglen_count == self.tracklength[i][j] and l == 0:
                            l += 1
                            seglen_count = 0
                            k -= 1
                        #elif count > 1 and seglen_count == self.tracklength[i][j] and l > 0:

                        if count > 3 and seglen_count < self.tracklength[i][j]:
                            k -= 2
                        if count > 3 and seglen_count == self.tracklength[i][j] and l > 0:
                            k += 3
                            continue
                        k += 1


                        """
                        if count > 1 and l < 1:
                            k -= 1
                            l += 1
                        elif count > 1 and l >= 1:
                            l = 0
                            k += count
                         """

                    """if k < self.num_segments:
                        continue
                    else:
                        break"""


    def scalarFlux(self):
        pass