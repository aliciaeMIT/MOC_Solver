import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class InitializeTracks(object):
    def __init__(self, num_azim, spacing, width, height, num_polar, radius, num_rings = 0, ring_radii = None):
        """
        This class generates tracks for method of characteristics, and their quadrature (azimuthal and polar).
        """
        self.num_azim2 = num_azim /2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.n_p = num_polar
        self.radius = radius
        self.phi = []
        self.nx = []
        self.ny = []
        self.ntot = []
        self.phi_eff = []
        self.phi_comp = []
        self.t_eff = []
        self.dx = []
        self.dy = []
        self.startpoint = [[] for _ in range(self.num_azim2)]
        self.endpoint = [[] for _ in range(self.num_azim2)]

        self.intersect1 = np.empty((num_azim,100), dtype = object)
        self.intersect2 = np.empty((num_azim,100), dtype = object)
        self.segstart = np.empty((num_azim * 100), dtype = object)
        self.segend = np.empty((num_azim * 100), dtype = object)
        self.seglength = np.empty((num_azim * 100), dtype = object)
        self.segangle = np.empty((num_azim * 100), dtype = object)
        self.segvolume = np.empty((num_azim * 100), dtype = object)
        self.segarea = np.empty((num_azim * 100), dtype = object)
        self.segmidpt = np.empty((num_azim * 100), dtype = object)
        self.segsource = np.empty((num_azim * 100), dtype = object)
        self.tracklengths = np.empty((num_azim, 100), dtype = object)
        self.boundids = np.empty((num_azim, 100), dtype = object)


    def getEnd(self, i, j):
        slope = math.tan(self.phi_eff[i])
        x0, y0 = self.startpoint[i][j]
        poss = []
        poss.append((0, y0 - slope * x0))
        poss.append((self.width, y0 + slope * (self.width - x0)))
        poss.append((x0 - y0 / slope, 0))
        poss.append((x0 - (y0 - self.height) / slope, self.height))

        for k in range(4):
            diffx = round(math.fabs(x0 - poss[k][0]),3)
            diffy = round(math.fabs(y0 - poss[k][1]),3)
            if (diffx == 0.0) and (diffy == 0.0):
                pass #checks for trivial case where it found the start point, does not save that as a solution
            elif (poss[k][0] >= 0.00 and poss[k][0] <= self.width) and (poss[k][1] >= 0.00 and poss[k][1] <=self.height):
                return poss[k]


    def makeTracks(self):
        self.tracks = []
        print "Getting ray entrance coordinates...\n"
        print "Getting ray exit coordinates...\n"
        for i in range(self.num_azim2):
            xin = np.zeros(self.ntot[i])
            yin = np.zeros(self.ntot[i])

            xin[:self.nx[i]] = self.dx[i]*(0.5 + np.arange(self.nx[i]))
            yin[:self.nx[i]] = 0
            yin[self.nx[i]:]= self.dy[i] * (0.5 + np.arange(self.ny[i]))

            if math.sin(self.phi_eff[i]) > 0 and math.cos(self.phi_eff[i]) > 0:
                #xin[self.nx:] = 0
                #redundant to set it to zero; initialized with zeros.
                pass
            elif math.sin(self.phi_eff[i]) > 0 and math.cos(self.phi_eff[i]) < 0:
                xin[self.nx[i]:] = self.width
            else:
                print "Error in makeTracks method"

            for j in range(int(self.ntot[i])):
                self.startpoint[i].append((xin[j], yin[j]))
                self.endpoint[i].append(self.getEnd(i,j))


    def plotTracks(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)
        ax1.add_patch(c)

        for i in range(self.num_azim2):
        #for i in [3, 7]:

            counter = int(self.ntot[i])

            for j in range(counter):
                try:
                    x1 = (self.startpoint[i][j][0])
                    x2 = (self.endpoint[i][j][0])
                    if x1 == x2:
                        print "Error! X values are equal for i = %d, j = %d" %(i,j)
                        print "x1 = %f \t x2 = %f" %(x1, x2)
                    y1 = self.startpoint[i][j][1]
                    y2 = self.endpoint[i][j][1]
                    if y1 == y2:
                        print "Error! y values are equal for i= %d, j = %f" %(i,j)
                        print "y1 = %f \t y2 = %f" %(y1, y2)
                    xvals = [x1, x2]
                    yvals = [y1,y2]
                except(TypeError):
                    print "i out of n_azim2", i, self.num_azim2
                    print "j out of counter", j, counter
                    raise

                plt.plot(xvals, yvals, 'k')

                if not (self.intersect1[i][j] == None):
                    xi1, yi1 = self.intersect1[i][j]
                    xi2, yi2 = self.intersect2[i][j]
                    plt.plot(xi1, yi1,'ro')
                    plt.plot(xi2, yi2, 'go')

        print "plotting tracks..."
        plt.show()

    def getTracks(self):
        print "\n------------------\nInput parameters:\n------------------"
        print "Num of azimuthal angles desired = %d" %(self.num_azim2 * 2)
        print "Track spacing = %.4f cm" %(self.spacing)
        print "Width of geometry = %.4f cm\nHeight of geometry = %.4f cm" %(self.width , self.height)
        print "\n"

        for i in range(self.num_azim2/2):
            if ((self.num_azim2 * 2)% 4 == 0):
                pass
            else:
                print "Error! Number of azimuthal angles must be a multiple of 4 for reflective geometry."
                break

            self.phi.append(math.pi / self.num_azim2 * (0.5 + i))
            self.nx.append(int(math.fabs((self.width / self.spacing) * math.sin(self.phi[i])) + 1))
            self.ny.append(int(math.fabs((self.height / self.spacing) * math.cos(self.phi[i]))))
            self.phi_eff.append(math.atan2((self.height * self.nx[i]), (self.width * self.ny[i])))
            self.phi_comp.append(math.pi - self.phi_eff[i])
            self.t_eff.append((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))
            self.dx.append(self.width/self.nx[i])
            self.dy.append(self.height / self.ny[i])


        for i in range(self.num_azim2/2):
            #complementary angle
            k = i
            self.phi_eff.append(self.phi_comp[k])
            self.nx.append(self.nx[k])
            self.ny.append(self.ny[k])
            self.dx.append(self.dx[k])
            self.dy.append(self.dy[k])
            self.phi.append(math.pi / self.num_azim2 * (0.5 + i))
            self.phi_comp.append(math.pi - self.phi_eff[i])
            self.t_eff.append((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))


        for i in range(self.num_azim2):
            self.ntot.append(self.nx[i] + self.ny[i])
            print "------------------\n Angle %d of %d \n------------------" % (i + 1, self.num_azim2)
            print "Phi = %f" % (math.degrees(self.phi[i]))
            print "nx = %f" % (self.nx[i])
            print "ny = %f" % (self.ny[i])
            print "ntot = %f" % (self.ntot[i])
            print "phi_eff = %f" % (math.degrees(self.phi_eff[i]))
            print "d_eff = %.3f cm" % (self.t_eff[i])
            print "dx = %.3f cm" % (self.dx[i])
            print "dy = %.3f cm" % (self.dy[i])
            print "phi_comp = %.3f" %(math.degrees(self.phi_comp[i]))
            print "\n"

    def getAngularQuadrature(self):
        """computation of azimuthal angle quadrature set, based on fraction of angular space of each angle.
        """
        omega_m_tot = 0
        self.omega_m = []
        for i in range(self.num_azim2):
            if (i == 0):
                self.omega_m.append((((self.phi[i + 1] - self.phi[i]) / 2) + self.phi[i]) / (2 * math.pi))
                omega_m_tot += self.omega_m[i]
            elif (i < (self.num_azim2 - 1)):
                self.omega_m.append((((self.phi[i + 1] - self.phi[i]) / 2) + ((self.phi[i] - self.phi[i - 1]) / 2)) / (2 * math.pi))
                omega_m_tot += self.omega_m[i]
            else:
                self.omega_m.append((2 * math.pi - self.phi[i] + (self.phi[i] - self.phi[i - 1]) / 2) / (2 * math.pi))
                omega_m_tot += self.omega_m[i]
        print "Calculating azimuthal weights...."
        print self.omega_m
        print "Total azimuthal weight sum: %f\n\n" %(omega_m_tot)
        return self.omega_m

    def getPolarWeight(self):
        """determines polar quadrature weights and angles using Tabuchi  Yamamoto (TY) set;
        number of polar divisions can be specified (2 or 3)
        returns ([polar weight per division], [sin-theta_polar per division])
        """
        polar_wt_total = 0
        print "Calculating polar weights..."
        if self.n_p == 2:
            self.omega_p = [0.212854 , 0.787146]
            self.sintheta_p = [0.363900 , 0.899900]
            polar_wt_total = self.omega_p[0] + self.omega_p[1]
            polar_angle_total = self.sintheta_p[0] + self.sintheta_p[1]
        elif self.n_p == 3:
            self.omega_p = [0.046233, 0.283619, 0.670148]
            self.sintheta_p = [0.166648, 0.537707, 0.932954]
            polar_wt_total = self.omega_p[0] + self.omega_p[1] + self.omega_p[2]
            polar_angle_total = self.sintheta_p[0] + self.sintheta_p[1] + self.sintheta_p[2]
        else:
            print "Error: must use 2 or 3 polar divisions for Tabuchi-Yamamoto polar quadrature."

        print "w_p:"
        print self.omega_p
        print "sin_theta_p:"
        print self.sintheta_p
        print "omega_p_total = %f" %(polar_wt_total)
        print "sintheta_p_total = %f\n\n" %(polar_angle_total)
        return self.omega_p, self.sintheta_p

    def findIntersection(self):
        self.num_segments=0 #increments index for storing segments
        print "Finding intersection points...\n\n"
        
        for i in range(self.num_azim2):
            for j in range(int(self.ntot[i])):
                cx0 = self.width/2
                cy0 = self.height/2
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
                    #distance from close to circle intersection point
                    dclose = math.sqrt(self.radius**2 - dist_center **2)

                    #first intersection point
                    fx = (close - dclose) * xproj + x0
                    fy = (close - dclose) * yproj + y0
                    self.intersect1[i][j] = (fx, fy)

                    #store first segment: from startpoint to intersect1
                    self.segmentStore(x0,fx, y0, fy, self.num_segments, i, j, 0)
                    #self.segmidpt[self.num_segments] = self.getMidpoint(x0, y0, fx, fy)
                    self.num_segments += 1 #increment to store next segment

                    #second intersection point
                    gx = (close + dclose) * xproj + x0
                    gy = (close + dclose) * yproj + y0
                    self.intersect2[i][j] = (gx, gy)

                    #store second segment: from intersect1 to intersect2
                    self.segmentStore(fx, gx, fy, gy, self.num_segments, i, j, 1)
                    #self.segmidpt[self.num_segments] = self.getMidpoint(fx, fy, gx, gy)
                    self.num_segments += 1 #increment to store next segment

                    #store third segment: from intersect2 to endpoint
                    self.segmentStore(gx, x1, gy, y1, self.num_segments, i, j, 0)
                    #self.segmidpt[self.num_segments] = self.getMidpoint(gx, gy, x1, y1)
                    self.num_segments += 1

                elif dist_center == self.radius:
                    print "line is tangent\n"
                    #treat as a miss. store whole track as 1 segment.
                    #later could improve this by calculating the point where it hits, segmenting into 2 at that point

                    self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    #self.segmidpt[self.num_segments] = self.getMidpoint(x0,y0, x1, y1)
                    self.num_segments += 1

                    #point e is tangent to circle; brushes but does not enter.
                else:
                    self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    #self.segmidpt[self.num_segments] = self.getMidpoint(x0,y0, x1, y1)
                    self.num_segments += 1

    def lengthTwoPoints(self, x1, x2, y1, y2):
        length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 )
        return length

    def getMidpoint(self, x1, y1, x2, y2):
        xc = abs((x2 - x1) / 2)
        yc = abs((y2-y1)/2)
        return (xc, yc)

    def segmentStore(self, x1, x2, y1, y2, k, i, j, q):
        self.segstart[k] = (x1, y1)
        self.segend[k] = (x2, y2)
        self.seglength[k] = self.lengthTwoPoints(x1, x2, y1, y2)
        self.segangle[k] = (i,j)  #store angle index i so phi, omega_m[i] can be retrieved later
        self.segsource[k] = q     #1 for fuel region; 0 for moderator


    def plotSegments(self):

        fig1 = plt.figure()

        ax1 = fig1.add_subplot(111, aspect='equal')

        plt.axis([0, self.width, 0, self.height])

        c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)

        ax1.add_patch(c)

        if not (self.num_rings == 0):

            rings = np.zeros(self.num_rings)

            for k in range(len(self.ring_radii)):

                ring = patches.Circle((self.width/2, self.height/2), (self.ring_radii[k]), fill=False)

                ax1.add_patch(ring)



        for k in range(self.num_segments):


            x1 = (self.segstart[k][0])

            x2 = (self.segend[k][0])

            if x1 == x2:

                print "Error! X values are equal for i = %d, j = %d" %(i,j)

                print "x1 = %f \t x2 = %f" %(x1, x2)

            y1 = self.segstart[k][1]

            y2 = self.segend[k][1]

            if y1 == y2:

                print "Error! y values are equal for i= %d, j = %f" %(i,j)

                print "y1 = %f \t y2 = %f" %(y1, y2)

            xvals = [x1, x2]

            yvals = [y1,y2]


            plt.plot(xvals, yvals)


        print "plotting segments..."

        plt.show()

    def getFSRAreas(self):
        area = 0
        quadweight = 0
        for k in range(self.num_segments):  #loop over all segments
            i = self.segangle[k][0]
            quadweight =  self.omega_m[i]  * self.t_eff[i]
            area +=  quadweight * self.seglength[k]
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



        for p in range(self.n_p):    #loop over polar angles
            for k in range(self.num_segments):  #loop over all segments
                i = self.segangle[k][0]
                self.segvolume[k] = self.omega_m[i] * self.t_eff[i] * self.seglength[k] *  self.sintheta_p[p]        

                quadweight +=  self.omega_m[i]  * self.t_eff[i] * self.omega_p[p]
                area +=  self.omega_m[i] * self.t_eff[i] * self.seglength[k]
                volume += self.segvolume[k]
        
        #estimated_volume = (4/3) * math.pi * (self.width / math.sqrt(2)) ** 3
        #expected_area = self.width * self.height     #area of pincell
        #print "volume calculated = %f \nvolume expected = %f \n" %(volume, estimated_volume)
        #print "area calculated = %f \nArea expected = %f \n" %(area, expected_area)

    def findBoundaryID(self, koords):
        #finds which boundary a start/endpoint lies on. takes in a tuple of koordinates
        #only applies for boundary points: where x=0 or x=xmax and y=0 or y=ymax.
        #boundary 1: bottom (y=0)
        #boundary 2: left (x=0)
        #boundary 3: top (y=ymax)
        #boundary 4: right (x=xmax)
        x, y = koords
        boundary_id = 0 #will remain 0 if not on a boundary.

        x = round(x, 3)
        y = round(y, 3)
        xmax = round(self.width,3)
        ymax = round(self.height, 3)
        #print "Finding boundary index..."
        #first identify if the point given lies on a boundary or is inside the box
        if not(x == 0.0 or x == xmax or y == 0.0 or y == ymax):
            print "Coordinate is not on boundary"
        else:

            if y == 0.0: #bottom
                boundary_id = 1
            elif x == 0.0: #left
                boundary_id = 2
            elif y == ymax: #top
                boundary_id = 3
            elif x == xmax:
                boundary_id = 4
            else:
                print "Error: boundary ID could not be determined! Check rounding/truncation of coordinates?"

        #print "Boundary ID for point (%.4f, %.4f): \t %d" %(x,y,boundary_id)
        return boundary_id


    def plotFluxPasses(self):

        fig1 = plt.figure()

        ax1 = fig1.add_subplot(111, aspect='equal')

        plt.axis([0, self.width, 0, self.height])

        #c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)

        #ax1.add_patch(c)

        tempvar = 1
        i = tempvar
        jmax = int(self.ntot[i])
        step = int(self.nx[i])
        loop1 =0
        jstart = 0
        for j in range(jstart, jmax, step):

            if loop1% 2 == 0:
                i = tempvar
            else:
                i = self.num_azim2 - i - 1


            x1 = (self.startpoint[i][j][0])

            x2 = (self.endpoint[i][j][0])

            y1 = self.startpoint[i][j][1]

            y2 = self.endpoint[i][j][1]

            xvals = [x1, x2]

            yvals = [y1,y2]


            plt.plot(xvals, yvals)
            loop1 += 1

        loop2 = loop1
        for j in range(int(self.ntot[tempvar] - jstart-1), -1*step, -1*step):
            if loop2%2 == 0:
                #i = self.num_azim2 - i - 1
                i = tempvar
            else:
                #i = tempvar
                i = self.num_azim2 - i - 1

            x1 = (self.startpoint[i][j][0])
            x2 = (self.endpoint[i][j][0])
            y1 = self.startpoint[i][j][1]
            y2 = self.endpoint[i][j][1]                             
            xvals = [x1, x2]
            yvals = [y1,y2]

            plt.plot(xvals, yvals)
            loop2 += 1


        plt.show()

    def plotScalarFlux(self, scalarflux):

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)
        ax1.add_patch(c)
        xvals = np.zeros(len(scalarflux))
        yvals = np.zeros(len(scalarflux))
        fluxes = np.zeros(len(scalarflux))
        koords = np.zeros(len(scalarflux))

        for k in range(self.num_segments):
            x1 = self.segmidpt[k][0]
            y1 = self.segmidpt[k][1]
            koords[k] = (x1, y1, k)
            xvals[k] = x1
            yvals[k] = y1
            fluxes[k] = scalarflux[k]
        
        #fluxes = np.reshape(scalarflux, [len(xvals), len(yvals)])
        #xvals, yvals = np.meshgrid(xvals,yvals)
        #fluxes = np.array(scalarflux)

        #plt.imshow(scalarflux)
        #plt.pcolormesh(xvals, yvals, fluxes.reshape(xvals.shape))
        heatmap, _, _ = np.histogram2d(xvals, yvals, weights=fluxes)
        plt.clf()
        plt.imshow(heatmap)
        plt.show()

class SingleTrack(object):
    def __init__(self, start_koords, end_koords, phi):
        """
        Class for creating a single track. Stores the incoming, outgoing coords,
        segments, incoming and outgoing track, angle
        """
        self.start_koords = start_koords
        self.end_koords = end_koords
        self.phi = phi
        self.track_in = None
        self.track_out = None
        self.segments = []
