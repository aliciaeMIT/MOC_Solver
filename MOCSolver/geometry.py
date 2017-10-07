import math
import numpy as np


class Geometry():
    def __init__(self, width, height, radius):
        self.width = width
        self.height = height
        self.radius = radius
        self.x0 = width / 2
        self.y0 = height / 2
        self.ring_radius = np.zeros(10) #store radii for rings to create in fuel



    def findIntersection(self, points_in, points_out):
        cx0 = self.width/2
        cy0 = self.height/2
        x0, y0 = points_in
        x1, y1 = points_out

        print "Finding intersection points..."

        raylen = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 )

        print "Ray length:"
        print raylen

        xproj = (x1 - x0) / raylen
        yproj = (y1 - y0) / raylen
        #consider putting a math.fabs() around these to prevent negative values if the math gets screwy?

        print "x, y projections:"
        print xproj
        print yproj

        close = xproj * (cx0 - x0) + yproj * (cy0 - y0)

        print "close to circle:"

        print close

        ex = xproj * close + x0
        ey = yproj * close + y0

        print "Line point closest to circle center:"

        print "(%f, \t %f)"%(ex, ey)

        dist_center = math.sqrt((ex - cx0) ** 2 + (ey - cy0) ** 2)

        print "distance to center from line:"

        print dist_center

        if dist_center < self.radius:
            #distance from close to circle intersection point
            dclose = math.sqrt(self.radius**2 - dist_center **2)

            #first intersection point
            fx = (close - dclose) * xproj + x0
            fy = (close - dclose) * yproj + y0

            #second intersection point
            gx = (close + dclose) * xproj + x0
            gy = (close + dclose) * yproj + y0
            print "Line intersects!"

            print "First point: (%.3f, %.3f)" %(fx, fy)
            print "Second point: (%.3f, %.3f" %(gx, gy)
            return (fx, fy), (gx, gy)

        elif dist_center == self.radius:
            print "line is tangent"
            return
            #point e is tangent to circle; brushes but does not enter.
        else:
            print "line does not intersect"
            return

    def createRings(self, num_rings):
        #fuel_edge = (x-x0) ** 2 + (y-y0) ** 2 - radius ** 2
        print "Creating rings..."

        ring_spacing = self.radius / (num_rings + 1)

        for i in range(0, num_rings):
            j = i + 1
            self.ring_radius[i] = self.radius - (j * ring_spacing)
            print "ring %d created: radius = %.3f cm" %(j, self.ring_radius[i])


        pass

    def divideSectors(self):
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