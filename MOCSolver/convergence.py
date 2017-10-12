import math
import numpy as np

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