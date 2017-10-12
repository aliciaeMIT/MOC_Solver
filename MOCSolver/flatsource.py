import math
import numpy as np


class FlatSourceApproximation(object):
    def __init__(self, q_fuel, q_mod, qweight):

        """
        This class generates the source for each flat source (and updates it, if not constant isotropic)
        """

        self.q_fuel = q_fuel
        self.q_mod = q_mod
        self.qweight = qweight
        self.qseg = np.zeros(1000)
        self.sigma_FSR = np.zeros(1000)


    def computeSource(self, num_seg, sigma_a, sigma_r):
        print "Determining source values for segment regions..."
        for k in range(num_seg):
            #self.qseg[k] = self.q_fuel * self.qweight[k]
            if self.qweight[k] == 0: #moderator
                self.sigma_FSR[k] = sigma_a
                self.qseg[k] = self.q_mod
            elif self.qweight[k] == 1:
                self.sigma_FSR[k] = sigma_r
                self.qseg[k] = self.q_fuel


    def checkVals(self):
        print self.qweight