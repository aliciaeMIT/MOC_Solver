import math
import numpy as np


class Geometry():
    def __init__(self, width, height, radius):
        self.width = width
        self.height = height
        self.radius = radius

    def createFuel(self, x, y):

        pass

    def divideRegions(self):
        pass

    def isFuelRegion(self):
        fxy = x ** 2 + y ** 2 - self.radius ** 2
        if fxy == 0:  # on circle edge
            return 1
        elif fxy > 0:  # in fuel rod
            # in fuel rod, have uniform isotropic source
            return 2
        elif fxy < 0:  # in moderator
            # in moderator, no source
            return 3