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

    def computeSource(self, num_seg):
        for k in range(num_seg):
            self.qseg[k] = self.q_fuel * self.qweight[k]


    def checkVals(self):
        print self.qweight