#Alicia M. Elliott, 22.212 Fall 2017
#Method of Characteristics solver
#2D pincell, fixed isotropic source in fuel

import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

class InitializeTracks(object):
    def __init__(self, num_azim, spacing, width, height, num_polar, radius, fsr, num_rings = 0, ring_radii = None):
        """
        This class generates tracks for method of characteristics, and their quadrature (azimuthal and polar).
        """
        self.num_azim2 = num_azim /2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.n_p = num_polar
        self.radius = radius
        self.fsr = fsr
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

        self.boundids = [[] for _ in range(self.num_azim2)]


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
        self.tracks = [[] for _ in range(self.num_azim2)] #self.tracks[i][j] to get a given track
        print("Getting ray entrance coordinates...\n")
        print("Getting ray exit coordinates...\n")

        # Create tracks in [0, pi/4]
        for i in range(self.num_azim2//2):
            for j in range(self.nx[i] + self.ny[i]):
                if j < self.nx[i]:
                    start = (self.width - (j+0.5) * self.dx[i], 0.0)
                else:
                    start = (0.0, (j-self.nx[i]+0.5) * self.dy[i])

                if j < self.ny[i]:
                    end = (self.width, (j+0.5)*self.dy[i])
                else:
                    end = (self.width - (j-self.ny[i]+0.5) * self.dx[i], self.height)

                self.startpoint[i].append(start)
                self.endpoint[i].append(end)
                thisTrack = SingleTrack(self.startpoint[i][j], self.endpoint[i][j], self.phi_eff[i])
                self.tracks[i].append(thisTrack)

        # Create tracks in [pi/4, pi/2]
        for i in range(self.num_azim2//2):
            ii = self.num_azim2//2 - i - 1
            for j in range(self.nx[ii] + self.ny[ii]):
                if j < self.nx[ii]:
                    start = ((j+0.5) * self.dx[ii], 0.0)
                else:
                    start = (self.width, (j-self.nx[ii]+0.5) * self.dy[ii])

                if j < self.ny[ii]:
                    end = (0.0, (j+0.5)*self.dy[ii])
                else:
                    end = ((j-self.ny[ii]+0.5) * self.dx[ii], self.height)

                self.startpoint[self.num_azim2//2+i].append(start)
                self.endpoint[self.num_azim2//2+i].append(end)
                thisTrack = SingleTrack(self.startpoint[self.num_azim2//2+i][j],
                                        self.endpoint[self.num_azim2//2+i][j],
                                        self.phi_eff[self.num_azim2//2+i])
                self.tracks[self.num_azim2//2+i].append(thisTrack)

    def plotTracks(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)
        ax1.add_patch(c)

        for i in range(self.num_azim2):
        #for i in [0, self.num_azim2-1]: #for debugging, to plot complementary angles (tracks should be cyclic)
            counter = int(self.ntot[i])
            for j in range(counter):
                try:
                    x1 = self.tracks[i][j].start_coords[0]
                    x2 = self.tracks[i][j].end_coords[0]

                    if x1 == x2:
                        print "Error! X values are equal for i = %d, j = %d" %(i,j)
                        print "x1 = %f \t x2 = %f" %(x1, x2)

                    y1 = self.tracks[i][j].start_coords[1]
                    y2 = self.tracks[i][j].end_coords[1]

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

    def getTrackParams(self):
        print "\n------------------\nInput parameters:\n------------------"
        print "Num of azimuthal angles desired = %d" %(self.num_azim2 * 2)
        print "Track spacing = %.4f cm" %(self.spacing)
        print "Width of geometry = %.4f cm\nHeight of geometry = %.4f cm" %(self.width , self.height)
        print "\n"

        for i in range(self.num_azim2/2):
            if ((self.num_azim2 * 2)% 4 == 0):
                pass
            else:
                print "Error: Number of azimuthal angles must be a multiple of 4 for reflective geometry."
                break

            self.phi.append(math.pi / self.num_azim2 * (0.5 + i))
            self.nx.append(int(np.ceil((self.width / self.spacing) * math.sin(self.phi[i]))))
            self.ny.append(int(np.ceil((self.height / self.spacing) * math.cos(self.phi[i]))))
            self.phi_eff.append(math.atan2((self.height * self.nx[i]), (self.width * self.ny[i])))
            self.phi_comp.append(math.pi - self.phi_eff[i])
            self.t_eff.append((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))
            self.dx.append(self.width/self.nx[i])
            self.dy.append(self.height / self.ny[i])


        for i in reversed(range(self.num_azim2//2)):
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
        """Calculation of azimuthal angle quadrature set, based on fraction of angular space of each angle.
        """
        omega_m_tot = 0
        self.omega_m = []
        for i in range(self.num_azim2):
            if (i == 0):
                self.omega_m.append(2*(((self.phi_eff[i + 1] - self.phi_eff[i]) / 2) + self.phi_eff[i]) / (2 * math.pi))
                omega_m_tot += self.omega_m[i]
            elif (i < (self.num_azim2 - 1)):
                self.omega_m.append(2*(((self.phi_eff[i + 1] - self.phi_eff[i]) / 2) + ((self.phi_eff[i] - self.phi_eff[i - 1]) / 2)) / (2 * math.pi))
                omega_m_tot += self.omega_m[i]
            else:
                self.omega_m.append(2*(math.pi - self.phi_eff[i] + (self.phi_eff[i] - self.phi_eff[i - 1]) / 2) / (2 * math.pi))
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
        self.intersect1 = [[] for _ in range(self.num_azim2)]
        self.intersect2 = [[] for _ in range(self.num_azim2)]
        print "Finding intersection points...\n\n"
        
        for i in range(self.num_azim2):
            for j in range(int(self.ntot[i])):
                cx0 = self.width/2
                cy0 = self.height/2
                x0, y0 = self.startpoint[i][j]
                x1, y1 = self.endpoint[i][j]
                track = self.tracks[i][j] #reference to object that stores this track
                raylen = self.lengthTwoPoints(x0, x1, y0, y1)

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
                    self.intersect1[i].append((fx, fy))

                    #store first segment: from startpoint to intersect1
                    s = self.segmentStore(x0,fx, y0, fy, self.num_segments, i, j, 0)
                    track.segments.append(s)

                    #second intersection point
                    gx = (close + dclose) * xproj + x0
                    gy = (close + dclose) * yproj + y0
                    self.intersect2[i].append((gx, gy))

                    #store second segment: from intersect1 to intersect2
                    s = self.segmentStore(fx, gx, fy, gy, self.num_segments, i, j, 1)
                    track.segments.append(s)

                    #store third segment: from intersect2 to endpoint
                    s = self.segmentStore(gx, x1, gy, y1, self.num_segments, i, j, 0)
                    track.segments.append(s)

                elif dist_center == self.radius:
                    print "line is tangent\n"
                    self.intersect1[i].append(None)
                    self.intersect2[i].append(None)
                    #treat as a miss. store whole track as 1 segment.
                    #later could improve this by calculating the point where it hits, segmenting into 2 at that point

                    s = self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    track.segments.append(s)

                    #point e is tangent to circle; brushes but does not enter.
                else:
                    self.intersect1[i].append(None)
                    self.intersect2[i].append(None)
                    s = self.segmentStore(x0, x1, y0, y1, self.num_segments, i, j, 0)
                    track.segments.append(s)

    def lengthTwoPoints(self, x1, x2, y1, y2):
        """finds distance between 2 points, returns the length"""
        length = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 )
        return length

    def segmentStore(self, x1, x2, y1, y2, k, i, j, q):
        start = (x1, y1)
        end = (x2,y2)
        newSeg = SingleSegment(start, end, q, self.lengthTwoPoints(x1, x2, y1, y2))
        return newSeg

    def reflectRays(self):

        print("Linking tracks...")
        for i in range(self.num_azim2//2):
            ii = self.num_azim2 - i - 1
            nx = self.nx[i]
            ny = self.ny[i]
            for j in range(nx + ny):
                track = self.tracks[i][j]

                # Set track in fwd direction
                if j < ny:
                    next_track = self.tracks[ii][nx+j]
                    track.track_in = next_track
                    track.refl_in = 1
                    next_track.track_out = track
                    next_track.refl_out = 0
                else:
                    next_track = self.tracks[ii][nx+2*ny-j-1]
                    track.track_in = next_track
                    track.refl_in = 0
                    next_track.track_in = track
                    next_track.refl_in = 0

                # Set track in bwd direction
                if j < nx:
                    next_track = self.tracks[ii][nx-j-1]
                    track.track_out = next_track
                    track.refl_out = 1
                    next_track.track_out = track
                    next_track.refl_out = 1
                else:
                    next_track = self.tracks[ii][j-nx]
                    track.track_out = next_track
                    track.refl_out = 0
                    next_track.track_in = track
                    next_track.refl_in = 1

    def plotSegments(self):
        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])

        #c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)
        #ax1.add_patch(c)

        for i in range(self.num_azim2):
            for track in self.tracks[i]:
                for s in track.segments:
                    x1, y1 = s.start_coords
                    x2, y2 = s.end_coords

                    if x1 == x2:
                        print "Error! X values are equal for i = %d, j = %d" %(i,j)
                        print "x1 = %f \t x2 = %f" %(x1, x2)

                    if y1 == y2:
                        print "Error! y values are equal for i= %d, j = %f" %(i,j)
                        print "y1 = %f \t y2 = %f" %(y1, y2)

                    xvals = [x1, x2]
                    yvals = [y1,y2]

                    if s.region == 0:
                        plt.plot(xvals, yvals, 'k')
                    elif s.region == 1:
                        plt.plot(xvals, yvals, 'r')
                    else:
                        print "Error: segment region not set"
        print "plotting segments..."
        plt.show()

    def getFSRVolumes(self, fuel, mod):
        """
            this is for computing FSR area/volumes and quadrature weights
            to get area only, set p range to 1. should get the area of the pincell out.
            To get total FSR volumes: summing over polar angles and all segments
        """
        print "Calculating FSR volumes..."
        
        area = 0
        quadweight = 0

        for p in range(self.n_p):                   #loop over polar angles
            for i in range(self.num_azim2):         #loop over all angles
                for track in self.tracks[i]:        #loop over all tracks
                    for s in track.segments:        #loop over all segments
                        s.volume = self.omega_m[i] * self.t_eff[i] * s.length *  self.sintheta_p[p]
                        quadweight = self.omega_m[i]  * self.t_eff[i] * self.omega_p[p]
                        s.area = quadweight * s.length
                        area += s.area
                        if s.region == 0:
                            mod.area += s.area
                        elif s.region == 1:
                            fuel.area += s.area

        est_area = self.width * self.height     #area of pincell
        est_area_fuel = math.pi * self.radius ** 2
        est_area_mod = est_area - est_area_fuel
        tot_area = fuel.area + mod.area

        print "fuel area calculated = %f \nfuel area expected = %f \n" %(fuel.area, est_area_fuel)
        print "mod area calculated = %f \nmod area expected = %f \n" %(mod.area, est_area_mod)
        print "pincell area calculated = %f \npincell area expected = %f\n" %(tot_area, est_area)

        print "Correcting track lengths...\n"

        corr_fuel = est_area_fuel / fuel.area
        corr_mod = est_area_mod / mod.area
        print "fuel track area correction factor: %f \nmod track area correction factor: %f\n" %(corr_fuel, corr_mod)

        for i in range(self.num_azim2):  # loop over all angles
            for track in self.tracks[i]:  # loop over all tracks
                for s in track.segments:  # loop over all segments
                    if s.region == 0:
                        s.length *= corr_mod

                    elif s.region == 1:
                        s.length *= corr_fuel

        fuel.area = est_area_fuel
        mod.area = est_area_mod

    def findBoundaryID(self, coords):
        #finds which boundary a start/endpoint lies on. takes in a tuple of coordinates
        #only applies for boundary points: where x=0 or x=xmax and y=0 or y=ymax.
        #boundary 1: bottom (y=0)
        #boundary 2: left (x=0)
        #boundary 3: top (y=ymax)
        #boundary 4: right (x=xmax)
        x, y = coords
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

    def plotTrackLinking(self):

        plt.figure(figsize=(12,9))
        fwd = True
        for i in range(1,13):
            plt.subplot(3,4,i)
            track = self.tracks[0][0]
            fwd = True
            for j in range(i):
                plt.plot([track.start_coords[0], track.end_coords[0]],
                         [track.start_coords[1], track.end_coords[1]], 'k')

                # Get the next track and bool indicating direction on next track
                if fwd:
                    fwd = track.refl_in
                    track = track.track_in
                else:
                    fwd = track.refl_out
                    track = track.track_out

            plt.title('num tracks = {}'.format(i))
            plt.xlim([0, self.width])
            plt.ylim([0, self.height])
        #plt.tight_layout()
        plt.show()
        #plt.savefig('connecting_tracks.png')
        #plt.close()

    def getTrackLinkCoords(self):
        for i in range(self.num_azim2):
            for track in self.tracks[i]:
                start = self.findBoundaryID(track.start_coords)
                end = self.findBoundaryID(track.end_coords)
                self.boundids.append((start, end))

        for i in range(self.num_azim2):
            j=0
            for track in self.tracks[i]:
                j+=1
                startvalx, startvaly = track.start_coords
                startvalx = round(startvalx, 3)
                startvaly = round(startvaly, 3)
                startval = (startvalx, startvaly)
                for g in range(self.num_azim2):
                    h=0
                    for track1 in self.tracks[g]:
                        h+=1
                        endvalx, endvaly = track1.end_coords
                        startx2, starty2 = track1.start_coords

                        endvalx = round(endvalx, 3)
                        endvaly = round(endvaly, 3)
                        endval = (endvalx, endvaly)

                        startx2 = round(startx2, 3)
                        starty2 = round(starty2, 3)
                        start2 = (startx2, starty2)

                        if startval == endval or startval == start2 and not (track == track1):
                            print "coordinates match! startval[%d][%d] == endval[%d][%d]" % (i, j, g, h)

    def plotScalarFlux(self, flux_fuel, flux_mod):

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111, aspect='equal')
        plt.axis([0, self.width, 0, self.height])
        c = patches.Circle((self.width/2, self.height/2), self.radius, color='b', fill=True)
        ax1.add_patch(c)

        xvals = np.linspace(0, self.width, 100)
        yvals = np.zeros(len(xvals))
        for i in range(len(xvals)):
            dist = self.width/2 - self.radius
            if xvals[i] < dist or xvals[i] > (dist + 2 * self.radius):
                yvals[i] = flux_mod
            else:
                yvals[i] = flux_fuel

        plt.plot(xvals, yvals, 'm--')
        plt.show()

        #fluxes = np.reshape(scalarflux, [len(xvals), len(yvals)])
        #xvals, yvals = np.meshgrid(xvals,yvals)
        #fluxes = np.array(scalarflux)

        #plt.imshow(scalarflux)
        #plt.pcolormesh(xvals, yvals, fluxes.reshape(xvals.shape))
        #heatmap, _, _ = np.histogram2d(xvals, yvals, weights=fluxes)
        #plt.clf()
        #plt.imshow(heatmap)
        #plt.show()


class SingleTrack(object):
    def __init__(self, start_coords, end_coords, phi):
        """
        Class for creating a single track. Stores the incoming, outgoing coords,
        segments, incoming and outgoing track, angle
        """
        self.start_coords = start_coords
        self.end_coords = end_coords
        self.phi = phi
        self.track_in = None
        self.track_out = None
        self.segments = []
        self.flux_in = np.zeros((2,3))


class SingleSegment(object):
    def __init__(self, start_coords, end_coords, region, length):

        self.start_coords = start_coords
        self.end_coords = end_coords
        self.region = region
        self.length = length
        self.exponential = []