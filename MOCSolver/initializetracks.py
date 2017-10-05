import math
import numpy as np

class initializeTracks(object):
    def __init__(self, num_azim, spacing, width, height):
        self.num_azim2 = num_azim/2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.phi = []
        self.nx = []
        self.ny = []
        self.phi_eff = []
        self.phi_comp = []
        self.t_eff = []
        self.dx = []
        self.dy = []
        self.startpoint = np.empty((16, 50), dtype=object)
        self.endpoint = np.empty((16, 50), dtype=object)


    def getStart(self, i):
        print "Getting ray entrance coordinates...\n"
        print "\tCalculating x-axis crossings...\n"

        for j in range(0, int(self.nx[i])):
            temp = (j * self.dx[i] + (self.dx[i]/2))
            self.startpoint[i][j] = (temp,0)
            print self.startpoint[i][j]

        print "\t\tCalculating y-axis crossings...\n"
        for j in range(0, int(self.ny[i])):
            temp = ((j * self.dy[i] + (self.dy[i] / 2)))
            self.startpoint[i][(self.nx[i] + j)] = (0,temp)
            print self.startpoint[i][(self.nx[i] + j)]

    def getEnd(self, i):
        print "Getting ray exit coordinates...\n"
        slope = math.tan(self.phi[i])
        print "\tCalculating x-axis crossings...\n"

        for j in range(0, int(self.nx[i])):
            temp = 0
            self.endpoint[i][j] = (temp,0)
            print self.endpoint[i][j]

        print "\t\tCalculating y-axis crossings...\n"
        for j in range(0, int(self.ny[i])):
            temp = 0
            self.endpoint[i][(self.nx[i] + j)] = (0,temp)
            print self.endpoint[i][(self.nx[i] + j)]

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
            self.phi.append(math.pi / self.num_azim2 * (0.5 + i))
            print "Phi = %f" %(math.degrees(self.phi[i]))
            self.nx.append(int(math.fabs((self.width / self.spacing) * math.sin(self.phi[i])) + 1))
            print "nx = %f" %(self.nx[i])
            self.ny.append(int(math.fabs((self.height / self.spacing) * math.cos(self.phi[i])) + 1))
            print "ny = %f" %(self.ny[i])
            self.phi_eff.append(math.atan((self.height * self.nx[i]) / (self.width * self.ny[i])))
            print "phi_eff = %f" % (math.degrees(self.phi_eff[i]))
            self.phi_comp.append(math.pi - self.phi_eff[i])
            print "phi_comp = %f" % (math.degrees(self.phi_comp[i]))
            self.t_eff.append((self.width / self.nx[i]) * math.sin(self.phi_eff[i]))
            print "t_eff = %.3f cm" % (self.t_eff[i])
            self.dx.append(self.width/self.nx[i])
            print "dx = %.3f cm" % (self.dx[i])
            self.dy.append(self.height/self.ny[i])
            print "dy = %.3f cm" % (self.dy[i])
            print "\n"

            self.getStart(i)