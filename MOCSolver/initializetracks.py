import math
import numpy as np
import matplotlib.pyplot as plt

class initializeTracks(object):
    def __init__(self, num_azim, spacing, width, height):
        self.num_azim2 = num_azim /2
        self.spacing = spacing
        self.width = width
        self.height = height
        self.phi = np.zeros(50)
        self.nx = np.zeros(50)
        self.ny = np.zeros(50)
        self.ntot = np.zeros(50)
        self.phi_eff = np.zeros(50)
        self.phi_comp = np.zeros(50)
        self.t_eff = np.zeros(50)
        self.dx = np.zeros(50)
        self.dy = np.zeros(50)
        self.deff = np.zeros(50)
        self.startpoint = np.empty((64, 50), dtype=object)
        self.xvals = np.zeros((64,40))
        self.yvals = np.zeros((64,40))
        self.endpoint = np.empty((64, 50), dtype=object)
        self.poss = np.empty(4, dtype = object)


    def getStart(self):
        print "Getting ray entrance coordinates...\n"
        print "\tCalculating nx rays...\n"
        for i in range(0,self.num_azim2):
            for j in range(0, int(self.nx[i])):
               if (i < self.num_azim2 / 2):
                   temp = self.dx[i] * (self.nx[i] - j - 0.5)
                   self.startpoint[i][j] = (temp, 0)
               else:
                   temp = (j * self.dx[i] + (self.dx[i] / 2))
                   self.startpoint[i][j] = (temp, 0)


            print "\t\tCalculating ny rays...\n"
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
            slope = math.tan(self.phi[i]) #each ray has same slope for a given angle. only compute once.
            #print slope
            print "\tCalculating crossings...\n"
            counter = int(self.nx[i] + self.ny[i])
         #   i=1
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
                    if (diffx == 0.0) and (diffy == 0.0):#self.poss[k] == (x0, y0):
                        pass
                    elif (self.poss[k][0] >= 0.00 and self.poss[k][0] <= self.width) and (self.poss[k][1] >= 0.00 and self.poss[k][1] <=self.height):
                        self.endpoint[i][j] = self.poss[k]
                        #break
                        #continue
                    else:
                        pass
                        #print "endpoint ", self.endpoint[i][j]     , i, j


    def plotSingleAngle(self):
        for i in range(0, self.num_azim2+1):
            counter = int(self.ntot[i])

            for z in range(counter):
                try:
                    x1 = (self.startpoint[i][z][0])
                    x2 = (self.endpoint[i][z][0])
                    if x1 == x2:
                        print "Error! X values are equal for i = %d, j = %d" %(i,z)
                        print "x1 = %f \t x2 = %f" %(x1, x2)
                    y1 = self.startpoint[i][z][1]
                    y2 = self.endpoint[i][z][1]
                    if y1 == y2:
                        print "Error! y values are equal for i= %d, j = %f" %(i,z)
                        print "y1 = %f \t y2 = %f" %(y1, y2)
                    xvals = [x1, x2]
                    yvals = [y1,y2]
                except(TypeError):
                    print "i out of n_azim2", i, self.num_azim2
                    print "z out of counter", z, counter
                    raise
                #plt.plot(self.startpoint[i][j],self.endpoint[i][j])
                plt.plot(xvals, yvals)
                #plt.plot(self.xvals[i][j], self.yvals[i][j])
                plt.axis([0, self.width, 0, self.height])
        print "plotting angle %d ray %d" %(i, z)
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
        self.plotSingleAngle()