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
        for i in range(len(vec_n)):
            sum1 += ((vec_n[i] - vec_n1[i])/vec_n[i]) ** 2
        I = len(vec_n)
        l2 = math.sqrt((1 / I) * sum1)
        if l2 < epsilon:
            return True
        else:
            return False