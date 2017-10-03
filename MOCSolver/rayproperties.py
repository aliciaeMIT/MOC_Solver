import math
import numpy as np
import matplotlib.pyplot as plt


class RayProperties():
    def __init__(self, phi_eff, ray_num, spacing_eff, cellwidth, cellheight):
        self.phi_eff = phi_eff
        self.ray_num = ray_num
        self.spacing_eff = spacing_eff
        self.width = cellwidth
        self.height = cellheight


    def slope(self):
        return (math.tan(self.phi_eff))

    def b(self, x_0, y_0, slope):
        if y_0 == 0:
            return (-1 * (slope * (self.ray_num) * self.spacing_eff))
        else:
            return (y_0 - slope * x_0)

    def entrance_coords(self):
        x = (self.ray_num+1) * self.spacing_eff
        return x, 0.0

    def find_exit_coords(self, slope, b, trace_num):
        if trace_num == 1:
            print "Ray trace 1"
            y = slope * self.width + b
            if (y <= self.height):
                return self.width, y
            else:
                x = (self.height - b) / slope
                if (x <= self.width):
                    return x, self.height
                else:
                    print "Error: check ray parameters. Cannot determine exit coordinates on trace 1."
        elif trace_num == 2:
            print "Ray trace 2"
            slope1 = -1 * (1 /slope)
            x = (self.height - b) / slope1
            if (x >= 0):
                return x, self.height
            else:
                y = slope1 * 0 + b
                if y<= self.height:
                    return y, 0
                else:
                    print "Error: check ray parameters. Cannot determine exit coordinates on trace 2."
        elif trace_num == 3:
            return 0

    def plot_ray(self, entrance, exit):
        x = [entrance[0], exit[0]]
        y = [entrance[1], exit[1]]
        print "Plotting ray..."
        plt.plot(x, y)
        plt.axis([0, self.width, 0, self.height])
        plt.show()
        return 0

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