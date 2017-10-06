import math
import numpy as np
import matplotlib.pyplot as plt

class InitializeTracks(object):
    def __init__(self, num_azim, spacing, width, height, num_polar):
        """
        This class generates tracks for method of characteristics, and their quadrature (azimuthal and polar).
        """
        self.num_azim2 = num_azim /2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.n_p = num_polar
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
        self.poss = np.empty(4, dtype = object)
        self.omega_m = np.zeros(self.num_azim2)
        self.omega_p = np.zeros(self.n_p)
        self.sintheta_p = np.zeros(self.n_p)

    def getStart(self):
        print "Getting  entrance coordinates...\n"

        for i in range(0,self.num_azim2):
            for j in range(0, int(self.nx[i])):
               if (i < self.num_azim2 / 2):
                   temp = self.dx[i] * (self.nx[i] - j - 0.5)
                   self.startpoint[i][j] = (temp, 0)
               else:
                   temp = (j * self.dx[i] + (self.dx[i] / 2))
                   self.startpoint[i][j] = (temp, 0)


            for j in range(0, int(self.ny[i])):
                k = self.nx[i] + j
                temp = ((j * self.dy[i] + (self.dy[i] / 2)))
                if (i < self.num_azim2 / 2):
                    self.startpoint[i][k] = (0.0,temp)

                else:
                    self.startpoint[i][k] = (self.width, temp)


    def getEnd(self):
        for i in range(0, self.num_azim2):
            print "Getting ray exit coordinates...\n"
            slope = math.tan(self.phi[i])

            counter = int(self.nx[i] + self.ny[i])

            for j in range(0, counter):
                x0, y0 = self.startpoint[i][j]
                b = y0 - slope * x0
                """
                if math.floor(b) == 0:
                    print "b = 0 for i = %d, j = %d " %(i,j)
                    print "b, %f" %(b)
                    print "x0, %f" %(x0)
                    print "y0, %f" %(y0)
                    print "slope, %f" %(slope)
                """
                y = slope * 0 + b
                y1 = slope * self.width + b

                x = (- b) / slope
                x1 = (self.height - b) / slope
                self.poss[0] = (0, y)
                self.poss[1] = (self.width, y1)
                self.poss[2] = (x, 0)
                self.poss[3] = (x1, self.height)

                for k in range(0,4):
                    diffx = round(math.fabs(x0 - self.poss[k][0]),3)
                    diffy = round(math.fabs(y0 - self.poss[k][1]),3)
                    if (diffx == 0.0) and (diffy == 0.0):
                        pass
                    elif (self.poss[k][0] >= 0.00 and self.poss[k][0] <= self.width) and (self.poss[k][1] >= 0.00 and self.poss[k][1] <=self.height):
                        self.endpoint[i][j] = self.poss[k]
                        #break
                        #continue
                    else:
                        pass


    def plotTracks(self):
        for i in range(0, self.num_azim2):
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
                plt.axis([0, self.width, 0, self.height])
        print "plotting..."
        plt.show()

    def getTracks(self):
        print "\n------------------\nInput parameters:\n------------------"
        print "Num of azimuthal angles desired = %d" %(self.num_azim2 * 2)
        print "Track spacing = %.4f cm" %(self.spacing)
        print "Width of geometry = %.4f cm\nHeight of geometry = %.4f cm" %(self.width , self.height)
        print "\n"

        for i in range(0, self.num_azim2):
            if ((self.num_azim2 * 2)% 4 == 0):
                pass
            else:
                print "Error! Number of azimuthal angles must be a multiple of 4."
                break

            print "------------------\n Angle %d of %d \n------------------" %(i+1, self.num_azim2)
            self.phi[i] = math.pi / self.num_azim2 * (0.5 + i)
            print "Phi = %f" %(math.degrees(self.phi[i]))
            self.nx[i] = int(math.fabs((self.width / self.spacing) * math.sin(self.phi[i])) + 1)
            print "nx = %f" %(self.nx[i])
            self.ny[i] = int(math.fabs((self.height / self.spacing) * math.cos(self.phi[i])) + 1)
            print "ny = %f" %(self.ny[i])
            self.ntot[i] = (self.nx[i] + self.ny[i])
            print "ntot = %f" %(self.ntot[i])
            self.phi_eff[i] = (math.atan((self.height * self.nx[i]) / (self.width * self.ny[i])))
            #self.phi_eff[i] = (math.atan2((self.height * self.nx[i]) , (self.width * self.ny[i])))
            print "phi_eff = %f" % (math.degrees(self.phi_eff[i]))
            self.phi_comp[i] = (math.pi - self.phi_eff[i])
            print "phi_comp = %f" % (math.degrees(self.phi_comp[i]))
            self.t_eff[i] = ((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))
            print "t_eff = %.3f cm" % (self.t_eff[i])
            self.dx[i] =(self.width/self.nx[i])
            print "dx = %.3f cm" % (self.dx[i])
            self.dy[i] = (self.height/self.ny[i])
            print "dy = %.3f cm" % (self.dy[i])
            self.deff[i] = (self.dx[i] * math.sin(self.phi_eff[i]))
            print "d_eff = %.3f cm" %(self.deff[i])
            print "\n"

            #complementary angle
            self.nx[self.num_azim2 - i -1] = self.nx[i]
            self.ny[self.num_azim2 - i -1] = self.ny[i]
            self.ntot[self.num_azim2 - i - 1] = self.ntot[i]
            self.dx[self.num_azim2 - i - 1] = self.dx[i]
            self.dy[self.num_azim2 - i - 1] = self.dy[i]
            self.deff[self.num_azim2 - i - 1] = self.deff[i]

            
        self.getStart()
        self.getEnd()
        self.plotTracks()
        self.getAngularQuadrature()
        self.getPolarWeight()

    def getAngularQuadrature(self):
        """computation of azimuthal angle quadrature set, based on fraction of angular space of each angle.
        """
        for i in range(self.num_azim2):
            if (i == 0):
                self.omega_m[i] = (((self.phi[i+1] - self.phi[i]) / 2) + self.phi[i]) / (2 * math.pi)
            elif (i < (self.num_azim2 - 1)):
                self.omega_m[i] = (((self.phi[i+1] - self.phi[i]) / 2) + ((self.phi[i] - self.phi[i-1])/ 2)) / (2 * math.pi)
            else:
                self.omega_m[i] = (2 * math.pi - self.phi[i] + (self.phi[i] - self.phi[i-1]) / 2)/ (2 * math.pi)
        print "Calculating azimuthal weights...."
        print self.omega_m
        return self.omega_m

    def getPolarWeight(self):
        """determines polar quadrature weights and angles using Tabuchi and Yamamoto (TY) set;
        number of polar divisions can be specified (2 or 3)
        returns ([polar weight per division], [sin-theta_polar per division])
        """
        print "Calculating polar weights..."
        if self.n_p == 2:
            self.omega_p = [0.212854 , 0.787146]
            self.sintheta_p = [0.363900 , 0.899900]

        elif self.n_p == 3:
            self.wp = [0.046233, 0.283619, 0.670148]
            self.sintheta_p = [0.166648, 0.537707, 0.932954]

        else:
            print "Error: must use 2 or 3 polar divisions for polar quadrature."

        print self.omega_p
        print self.sintheta_p
        return self.omega_p, self.sintheta_p