import math
import numpy as np
import matplotlib.pyplot as plt


class RayProperties():
    def __init__(self, phi, phi_eff, ray_num, spacing_eff, dx, dy, cellwidth, cellheight):
        self.phi = phi
        self.phi_eff = phi_eff
        self.ray_num = ray_num
        self.spacing_eff = spacing_eff
        self.dx = dx
        self.dy = dy
        self.width = cellwidth
        self.height = cellheight


    def slope(self, num_trace):
        #if num_trace % 2 == 1:
        return (math.tan(self.phi_eff))
        #elif num_trace % 2 == 0:
        #    return (math.tan(self.phi_eff))

    def b(self, x_0, y_0, slope, num_trace):
        if (num_trace % 2 == 1):
            b = (y_0 - slope * x_0)
            #b = (-1 * (slope * (self.ray_num) * self.spacing_eff))
            return b
        else:
            slope1 = -1 * (1 / slope)
            return (y_0 - slope1 * x_0)

    def entrance_coords(self):
        x = self.ray_num * self.dx +  (self.dx / 2)
        return x, 0.0

    def find_exit_coords(self, slope, b, trace_num):
        if (trace_num == 1) or (trace_num == 3):
            print "Ray trace %d" %(trace_num)
            if trace_num == 1:
                x = self.width
            elif trace_num == 3:
                x = 0
            y = slope * x + b
            if (y <= self.height):
                print "y less than height"
                return x, y
            else:
                x = (self.height - b) / slope
                if (x <= self.width):
                    print "x less than width"
                    return x, y
                else:
                    print "Error: check ray parameters. Cannot determine exit coordinates on trace 1."
        elif (trace_num == 2) or (trace_num == 4):
            print "Ray trace %d" % (trace_num)
            slope1 = -1 * (1 /slope)
            if trace_num == 2:
                y = self.height
            elif trace_num == 4:
                y = 0
            x = (y - b) / slope1
            if (x >= 0):
                print "x less than width, greater than zero"
                return x, y
            else:
                y = slope1 * x + b
                print slope1
                print b
                if y<= self.height:
                    print "y less than height, x less than 0"
                    return x, y
                else:
                    print "Error: check ray parameters. Cannot determine exit coordinates on trace 2."


    def plot_ray(self, entrance, exit, exit2, exit3, exit4):
        x1 = [entrance[0], exit[0]]
        y1 = [entrance[1], exit[1]]
        x2 = [exit[0], exit2[0]]
        y2 = [exit[1], exit2[1]]
        x3 = [exit2[0], exit3[0]]
        y3 = [exit2[1], exit3[1]]
        x4 = [exit3[0], exit4[0]]
        y4 = [exit3[1], exit4[1]]
        print "Plotting ray..."
        plt.plot(x1, y1)
        plt.plot(x2, y2)
        plt.plot(x3, y3)
        plt.plot(x4, y4)
        plt.axis([0, self.width, 0, self.height])
        plt.show()
        return 0

    def plot_all_tracks(self, nx, ray_coords_in, ray_coords_out):
        nx = int(nx)
        for j in range(0, nx):
            for num_trace in range(0, 4):
                print "plotting..."
                plt.plot(ray_coords_in[j, num_trace], ray_coords_out[j, num_trace])
                plt.axis([0, self.width, 0, self.height])
        plt.show()

    def print_coords(self, entrance, exit):
        print "Entrance coordinates:"
        print entrance
        print "Exit coordinates:"
        print exit
        print "\n"

    def xlength(self):
        return (self.width - (self.ray_num * self.spacing_eff))


    def ylength(self, x_length):
        return (x_length * math.tan(self.phi_eff))