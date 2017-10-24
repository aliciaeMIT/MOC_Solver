import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


class InitializeTracks(object):
    def __init__(self, num_azim, spacing, width, height, num_polar, radius, num_rings, ring_radii):
        """
        This class generates tracks for method of characteristics, and their quadrature (azimuthal and polar).

        Attributes:
            num_azim2       =  number of azimuthal angles in [0,pi]
            spacing         = desired track spacing
            width           = width of pincell
            height          = height of pincell

        """
        self.num_azim2 = num_azim / 2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.n_p = num_polar
        self.radius = radius
        self.ring_radii = ring_radii
        self.num_rings = num_rings

        self.phi = np.zeros(100)
        self.nx = np.zeros(100)
        self.ny = np.zeros(100)
        self.ntot = np.zeros(100)
        self.phi_eff = np.zeros(100)
        self.phi_comp = np.zeros(100)
        self.t_eff = np.zeros(100)
        self.dx = np.zeros(100)
        self.dy = np.zeros(100)
        self.deff = np.zeros(100)
        self.startpoint = np.empty((num_azim, 100), dtype=object)
        self.endpoint = np.empty((num_azim, 100), dtype=object)
        self.poss = np.empty(4, dtype=object)
        self.omega_m = np.zeros(self.num_azim2)
        self.omega_p = np.zeros(self.n_p)
        self.sintheta_p = np.zeros(self.n_p)
        self.intersect1 = np.empty((num_azim, 100), dtype=object)
        self.intersect2 = np.empty((num_azim, 100), dtype=object)
        self.segstart = np.empty((num_azim * 100), dtype=object)
        self.segend = np.empty((num_azim * 100), dtype=object)
        self.seglength = np.empty((num_azim * 100), dtype=object)
        self.segangle = np.empty((num_azim * 100), dtype=object)
        self.segvolume = np.empty((num_azim * 100), dtype=object)
        self.segarea = np.empty((num_azim * 100), dtype=object)
        self.segmidpt = np.empty((num_azim * 100), dtype=object)
        self.segsource = np.empty((num_azim * 100), dtype=object)
        self.tracklengths = np.empty((num_azim, 100), dtype=object)
        self.boundids = np.empty((num_azim, 100), dtype=object)

    def getStart(self):
        print "Getting ray entrance coordinates...\n"

        for i in range(0, self.num_azim2):
            for j in range(0, int(self.nx[i])):
                if (i >= self.num_azim2 / 2):
                    temp = self.dx[i] * (self.nx[i] - j - 0.5)
                    self.startpoint[i][j] = (temp, 0.0)
                else:
                    temp = (j * self.dx[i] + (self.dx[i] / 2))
                    self.startpoint[i][j] = (temp, 0.0)

            for j in range(0, int(self.ny[i])):
                k = self.nx[i] + j
                temp = j * self.dy[i] + (self.dy[i] / 2)
                if (i >= self.num_azim2 / 2):
                    self.startpoint[i][k] = (0.0, temp)
                else:
                    self.startpoint[i][k] = (self.width, temp)

    def getEnd(self):
        print "Getting ray exit coordinates...\n"
        for i in range(0, self.num_azim2):
            slope = math.tan(self.phi_eff[i])
            counter = int(self.nx[i] + self.ny[i])
            for j in range(0, counter):
                x0, y0 = self.startpoint[i][j]
                self.poss[0] = (0, y0 - slope * x0)
                self.poss[1] = (self.width, y0 + slope * (self.width - x0))
                self.poss[2] = (x0 - y0 / slope, 0)
                self.poss[3] = (x0 - (y0 - self.height) / slope, self.height)

                for k in range(0, 4):
                    diffx = round(math.fabs(x0 - self.poss[k][0]), 3)
                    diffy = round(math.fabs(y0 - self.poss[k][1]), 3)
                    if (diffx == 0.0) and (diffy == 0.0):
                        pass  # checks for trivial case where it found the start point, does not save that as a solution
                    elif (self.poss[k][0] >= 0.00 and self.poss[k][0] <= self.width) and (
                            self.poss[k][1] >= 0.00 and self.poss[k][1] <= self.height):
                        self.endpoint[i][j] = self.poss[k]


    def plotTracks(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width / 2, self.height / 2), self.radius, color='b', fill=True)
        ax1.add_patch(c)
 """
        if not (self.num_rings == 0):
            rings = np.zeros(self.num_rings)
            for k in range(len(self.ring_radii)):
                ring = patches.Circle((self.width / 2, self.height / 2), (self.ring_radii[k]), fill=False)
                ax1.add_patch(ring)
"""

        # for i in range(0, self.num_azim2):
        for i in [2, self.num_azim2 - 1 - 2]:
            counter = int(self.ntot[i])
            for j in range(counter):
                try:
                    x1 = (self.startpoint[i][j][0])
                    x2 = (self.endpoint[i][j][0])
                    if x1 == x2:
                        print "Error! X values are equal for i = %d, j = %d" % (i, j)
                        print "x1 = %f \t x2 = %f" % (x1, x2)
                    y1 = self.startpoint[i][j][1]
                    y2 = self.endpoint[i][j][1]
                    if y1 == y2:
                        print "Error! y values are equal for i= %d, j = %f" % (i, j)
                        print "y1 = %f \t y2 = %f" % (y1, y2)
                    xvals = [x1, x2]
                    yvals = [y1, y2]
                except(TypeError):
                    print "i out of n_azim2", i, self.num_azim2
                    print "j out of counter", j, counter
                    raise
                plt.plot(xvals, yvals, 'k')
                if not (self.intersect1[i][j] == None):
                    xi1, yi1 = self.intersect1[i][j]
                    xi2, yi2 = self.intersect2[i][j]
                    plt.plot(xi1, yi1, 'ro')
                    plt.plot(xi2, yi2, 'go')
        print "plotting tracks..."
        plt.show()

    def getTracks(self):
        print "\n------------------\nInput parameters:\n------------------"
        print "Num of azimuthal angles desired = %d" % (self.num_azim2 * 2)
        print "Track spacing = %.4f cm" % (self.spacing)
        print "Width of geometry = %.4f cm\nHeight of geometry = %.4f cm" % (self.width, self.height)
        print "\n"

        for i in range(0, self.num_azim2):
            if ((self.num_azim2 * 2) % 4 == 0):
                pass
            else:
                print "Error; Number of azimuthal angles must be a multiple of 4."
                break

            self.phi[i] = math.pi / self.num_azim2 * (0.5 + i)
            self.nx[i] = int(math.fabs((self.width / self.spacing) * math.sin(self.phi[i])) + 1)
            self.ny[i] = int(math.fabs((self.height / self.spacing) * math.cos(self.phi[i])))
            self.ntot[i] = (self.nx[i] + self.ny[i])
            self.phi_eff[i] = (math.atan((self.height * self.nx[i]) / (self.width * self.ny[i])))
            self.phi_comp[i] = (math.pi - self.phi_eff[i])  # self.phi_eff[i] + (math.pi/2)
            self.t_eff[i] = ((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))
            self.dx[i] = (self.width / self.nx[i])
            self.dy[i] = (self.height / self.ny[i])
            self.deff[i] = (self.dx[i] * math.sin(self.phi_eff[i]))
            # complementary angle
            comp_index = self.num_azim2 - i - 1
            self.phi_eff[comp_index] = self.phi_comp[i]
            self.nx[comp_index] = self.nx[i]
            self.ny[comp_index] = self.ny[i]
            self.ntot[comp_index] = self.ntot[i]
            self.dx[comp_index] = self.dx[i]
            self.dy[comp_index] = self.dy[i]

            print "------------------\n Angle %d of %d \n------------------" % (i + 1, self.num_azim2)
            print "Phi = %f" % (math.degrees(self.phi[i]))
            print "nx = %f" % (self.nx[i])
            print "ny = %f" % (self.ny[i])
            print "ntot = %f" % (self.ntot[i])
            print "phi_eff = %f" % (math.degrees(self.phi_eff[i]))
            print "t_eff = %.3f cm" % (self.t_eff[i])
            print "dx = %.3f cm" % (self.dx[i])
            print "dy = %.3f cm" % (self.dy[i])
            print "d_eff = %.3f cm" % (self.deff[i])
            print "\n"

    def findIntersection(self):
        self.num_segments = 0  # increments index for storing segments
        print "Finding intersection points...\n\n"

        for i in range(0, self.num_azim2):
            for j in range(int(self.ntot[i])):
                cx0 = self.width / 2
                cy0 = self.height / 2
                x0, y0 = self.startpoint[i][j]
                x1, y1 = self.endpoint[i][j]

                raylen = self.lengthTwoPoints(x0, x1, y0, y1)
                self.tracklengths[i][j] = raylen

                xproj = (x1 - x0) / raylen
                yproj = (y1 - y0) / raylen

                close = xproj * (cx0 - x0) + yproj * (cy0 - y0)

                ex = xproj * close + x0
                ey = yproj * close + y0

                dist_center = math.sqrt((ex - cx0) ** 2 + (ey - cy0) ** 2)

                if dist_center < self.radius:
                    # distance from close to circle intersection point
                    dclose = math.sqrt(self.radius ** 2 - dist_center ** 2)

                    # first intersection point
                    fx = (close - dclose) * xproj + x0
                    fy = (close - dclose) * yproj + y0
                    self.intersect1[i][j] = (fx, fy)

                    # store first segment: from startpoint to intersect1
                    self.segmentStore(x0, fx, y0, fy, self.num_segments, i, j, 0)
                    self.segmidpt[self.num_segments] = self.getMidpoint(x0, y0, fx, fy)
                    self.num_segments += 1  # increment to store next segment

                    # second intersection point
                    gx = (close + dclose) * xproj + x0
                    gy = (close + dclose) * yproj + y0
                    self.intersect2[i][j] = (gx, gy)

                    # store second segment: from intersect1 to intersect2
                    self.segmentStore(fx, gx, fy, gy, self.num_segments, i, j, 1)
                    self.segmidpt[self.num_segments] = self.getMidpoint(fx, fy, gx, gy)
                    self.num_segments += 1  # increment to store next segment

                    # store third segment: from intersect2 to endpoint
                    self.segmentStore(gx, x1, gy, y1, self.num_segments, i, j, 0)
                    self.segmidpt[self.num_segments] = self.getMidpoint(gx, gy, x1, y1)
                    self.num_segments += 1

                    """
                    print "Line intersects!"

                    print "First point: (%.3f, %.3f)" %(fx, fy)
                    print "Second point: (%.3f, %.3f\n" %(gx, gy)
                    """
                elif dist_center == self.radius:
                    print "line is tangent\n"
                    # treat as a miss. store whole track as 1 segment.
                    # later could improve this by calculating the point where it hits, segmenting into 2 at that point

                    self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    self.segmidpt[self.num_segments] = self.getMidpoint(x0, y0, x1, y1)
                    self.num_segments += 1

                    # point e is tangent to circle; brushes but does not enter.
                else:
                    self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    self.segmidpt[self.num_segments] = self.getMidpoint(x0, y0, x1, y1)
                    self.num_segments += 1
                    # print "line does not intersect"

    def lengthTwoPoints(self, x1, x2, y1, y2):
        length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)
        return length

    def getMidpoint(self, x1, y1, x2, y2):
        xc = abs((x2 - x1) / 2)
        yc = abs((y2 - y1) / 2)
        return (xc, yc)

    def segmentStore(self, x1, x2, y1, y2, k, i, j, q):
        self.segstart[k] = (x1, y1)
        self.segend[k] = (x2, y2)
        self.seglength[k] = self.lengthTwoPoints(x1, x2, y1, y2)
        self.segangle[k] = (i, j)  # store angle index i so phi, omega_m[i] can be retrieved later
        self.segsource[k] = q  # 1 for fuel region; 0 for moderator

    def plotSegments(self):

        fig1 = plt.figure()

        ax1 = fig1.add_subplot(111, aspect='equal')

        plt.axis([0, self.width, 0, self.height])

        c = patches.Circle((self.width / 2, self.height / 2), self.radius, color='b', fill=True)

        ax1.add_patch(c)

        if not (self.num_rings == 0):

            rings = np.zeros(self.num_rings)

            for k in range(len(self.ring_radii)):
                ring = patches.Circle((self.width / 2, self.height / 2), (self.ring_radii[k]), fill=False)

                ax1.add_patch(ring)

        for k in range(self.num_segments):

            x1 = (self.segstart[k][0])

            x2 = (self.segend[k][0])

            if x1 == x2:
                print "Error! X values are equal for i = %d, j = %d" % (i, j)

                print "x1 = %f \t x2 = %f" % (x1, x2)

            y1 = self.segstart[k][1]

            y2 = self.segend[k][1]

            if y1 == y2:
                print "Error! y values are equal for i= %d, j = %f" % (i, j)

                print "y1 = %f \t y2 = %f" % (y1, y2)

            xvals = [x1, x2]

            yvals = [y1, y2]

            plt.plot(xvals, yvals)

        print "plotting segments..."

        plt.show()

    def getFSRAreas(self):
        area = 0

        quadweight = 0
        for k in range(self.num_segments):  # loop over all segments
            i = self.segangle[k][0]
            quadweight = self.omega_m[i] * self.t_eff[i]
            area += quadweight * self.seglength[k]
            self.segarea[k] = quadweight * self.seglength[k]

    def getFSRVolumes(self):

        print "Calculating FSR volumes..."
        """
        this is for computing FSR area/volumes and quadrature weights
        to get area only, set p range to 1. should get the area of the pincell out.
        To get total FSR volumes: summing over polar angles and all segments
        """
        area = 0
        volume = 0
        quadweight = 0

        for p in range(self.n_p):  # loop over polar angles
            for k in range(self.num_segments):  # loop over all segments
                i = self.segangle[k][0]
                self.segvolume[k] = self.omega_m[i] * self.t_eff[i] * self.seglength[k] * self.sintheta_p[p]

                quadweight += self.omega_m[i] * self.t_eff[i] * self.omega_p[p]
                area += self.omega_m[i] * self.t_eff[i] * self.seglength[k]
                volume += self.segvolume[k]

                # estimated_volume = (4/3) * math.pi * (self.width / math.sqrt(2)) ** 3
                # expected_area = self.width * self.height     #area of pincell
                # print "volume calculated = %f \nvolume expected = %f \n" %(volume, estimated_volume)
                # print "area calculated = %f \nArea expected = %f \n" %(area, expected_area)

    def findBoundaryID(self, coords):
        # finds which boundary a start/endpoint lies on. takes in a tuple of coordinates
        # only applies for boundary points: where x=0 or x=xmax and y=0 or y=ymax.
        # boundary 1: bottom (y=0)
        # boundary 2: left (x=0)
        # boundary 3: top (y=ymax)
        # boundary 4: right (x=xmax)
        x, y = coords
        boundary_id = 0  # will remain 0 if not on a boundary.

        x = round(x, 3)
        y = round(y, 3)
        xmax = round(self.width, 3)
        ymax = round(self.height, 3)
        # print "Finding boundary index..."
        # first identify if the point given lies on a boundary or is inside the box
        if not (x == 0.0 or x == xmax or y == 0.0 or y == ymax):
            print "Coordinate is not on boundary"
        else:

            if y == 0.0:  # bottom
                boundary_id = 1
            elif x == 0.0:  # left
                boundary_id = 2
            elif y == ymax:  # top
                boundary_id = 3
            elif x == xmax:
                boundary_id = 4
            else:
                print "Error: boundary ID could not be determined! Check rounding/truncation of coordinates?"

        # print "Boundary ID for point (%.4f, %.4f): \t %d" %(x,y,boundary_id)
        return boundary_id

    def plotFluxPasses(self):

        fig1 = plt.figure()

        ax1 = fig1.add_subplot(111, aspect='equal')

        plt.axis([0, self.width, 0, self.height])

        # c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)

        # ax1.add_patch(c)

        tempvar = 1
        i = tempvar
        jmax = int(self.ntot[i])
        step = int(self.nx[i])
        loop1 = 0
        jstart = 0
        for j in range(jstart, jmax, step):

            if loop1 % 2 == 0:
                i = tempvar
            else:
                i = self.num_azim2 - i - 1

            x1 = (self.startpoint[i][j][0])

            x2 = (self.endpoint[i][j][0])

            y1 = self.startpoint[i][j][1]

            y2 = self.endpoint[i][j][1]

            xvals = [x1, x2]

            yvals = [y1, y2]

            plt.plot(xvals, yvals)
            loop1 += 1

        loop2 = loop1
        for j in range(int(self.ntot[tempvar] - jstart - 1), -1 * step, -1 * step):
            if loop2 % 2 == 0:
                # i = self.num_azim2 - i - 1
                i = tempvar
            else:
                # i = tempvar
                i = self.num_azim2 - i - 1

            x1 = (self.startpoint[i][j][0])
            x2 = (self.endpoint[i][j][0])
            y1 = self.startpoint[i][j][1]
            y2 = self.endpoint[i][j][1]
            xvals = [x1, x2]
            yvals = [y1, y2]

            plt.plot(xvals, yvals)
            loop2 += 1

        plt.show()

    def plotScalarFlux(self, scalarflux):

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width / 2, self.height / 2), self.radius, color='b', fill=True)
        ax1.add_patch(c)
        xvals = np.zeros(len(scalarflux))
        yvals = np.zeros(len(scalarflux))
        fluxes = np.zeros(len(scalarflux))
        coords = np.zeros(len(scalarflux))

        for k in range(self.num_segments):
            x1 = self.segmidpt[k][0]
            y1 = self.segmidpt[k][1]
            coords[k] = (x1, y1, k)
            xvals[k] = x1
            yvals[k] = y1
            fluxes[k] = scalarflux[k]

        # fluxes = np.reshape(scalarflux, [len(xvals), len(yvals)])
        # xvals, yvals = np.meshgrid(xvals,yvals)
        # fluxes = np.array(scalarflux)

        # plt.imshow(scalarflux)
        # plt.pcolormesh(xvals, yvals, fluxes.reshape(xvals.shape))
        heatmap, _, _ = np.histogram2d(xvals, yvals, weights=fluxes)
        plt.clf()
        plt.imshow(heatmap)
        plt.show()


    class Quadrature(object):
        """
        Polar quadrature set (Tabuchi-Yamamoto)

        Attributes:
            sintheta_p        = sin of polar angle
            omega_p           = polar weight
            omega_m           = azimuthal quadrature

        Methods:
            getPolarQuadrature      = sets polar quadrature from Tabuchi-Yamamoto
            getAngularQuadrature    = sets angular quadrature from input parameters
        """
        def __init__(self):
            self.sintheta_p = []
            self.omega_p = []

        def getAngularQuadrature(self, num_azim2, ):
            """computation of azimuthal angle quadrature set, based on fraction of angular space of each angle.
            """
            omega_m_tot = 0
            for i in range(num_azim2):
                """  this commented out block does not allow the sum of all quadrature weights to equal 2pi. 
                if (i == 0):
                    self.omega_m[i] = (((self.phi_eff[i+1] - self.phi_eff[i]) / 2) + self.phi_eff[i]) / (2 * math.pi)
                    omega_m_tot += self.omega_m[i]
                elif (i < (self.num_azim2-1)):
                    self.omega_m[i] = (((self.phi_eff[i+1] - self.phi_eff[i]) / 2) + ((self.phi_eff[i] - self.phi_eff[i-1])/ 2)) / (2 * math.pi )
                    omega_m_tot += self.omega_m[i]
                else:
                    self.omega_m[i] = (2 * math.pi - self.phi_eff[i] + (self.phi_eff[i] - self.phi_eff[i-1]) / 2)   / (2 * math.pi )
                    omega_m_tot += self.omega_m[i]
                """
                if (i == 0):
                    self.omega_m[i] = (((self.phi[i + 1] - self.phi[i]) / 2) + self.phi[i]) / (2 * math.pi)
                    omega_m_tot += self.omega_m[i]
                elif (i < (self.num_azim2 - 1)):
                    self.omega_m[i] = (((self.phi[i + 1] - self.phi[i]) / 2) + (
                    (self.phi[i] - self.phi[i - 1]) / 2)) / (
                                          2 * math.pi)
                    omega_m_tot += self.omega_m[i]
                else:
                    self.omega_m[i] = (2 * math.pi - self.phi[i] + (self.phi[i] - self.phi[i - 1]) / 2) / (
                    2 * math.pi)
                    omega_m_tot += self.omega_m[i]
                    # """
            print "Calculating azimuthal weights...."
            print self.omega_m
            print "Total azimuthal weight sum: %f\n\n" % (omega_m_tot)
            return self.omega_m

        def getPolarWeight(self):
            """determines polar quadrature weights and angles using Tabuchi and Yamamoto (TY) set;
            number of polar divisions can be specified (2 or 3)
            returns ([polar weight per division], [sin-theta_polar per division])
            """
            polar_wt_total = 0
            print "Calculating polar weights..."
            if self.n_p == 2:
                self.omega_p = [0.212854, 0.787146]
                self.sintheta_p = [0.363900, 0.899900]
                polar_wt_total = self.omega_p[0] + self.omega_p[1]
                polar_angle_total = self.sintheta_p[0] + self.sintheta_p[1]
            elif self.n_p == 3:
                self.omega_p = [0.046233, 0.283619, 0.670148]
                self.sintheta_p = [0.166648, 0.537707, 0.932954]
                polar_wt_total = self.omega_p[0] + self.omega_p[1] + self.omega_p[2]
                polar_angle_total = self.sintheta_p[0] + self.sintheta_p[1] + self.sintheta_p[2]
            else:
                print "Error: must use 2 or 3 polar divisions for polar quadrature."

            print "w_p:"
            print self.omega_p
            print "sin_theta_p:"
            print self.sintheta_p
            print "omega_p_total = %f" % (polar_wt_total)
            print "sintheta_p_total = %f\n\n" % (polar_angle_total)
            return self.omega_p, self.sintheta_p

    def ProblemSetup():