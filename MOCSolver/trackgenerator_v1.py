import math
import numpy as np

class TrackGeneration():
    def __init__(self, n_azim, spacing, width, height):
        self.n_azim = n_azim
        self.spacing = spacing
        self.width = width
        self.height = height



    def track_params(n_azim, spacing, width, height):

        # type: (int, float, float, float) -> None
        """This function outputs the effective azimuthal angle, track spacing, and azim angle quadrature weights
        """

        m = n_azim/4
        phi_eff = np.zeros(m)
        spacing_eff = np.zeros(m)

        if (n_azim % 4 == 0): #make sure m is a multiple of 4 before computing parameters
            print "Parameters input: \n number of azimuthal angles = %d \t track spacing = %f cm\n " \
                  "pincell width = %f cm \t pincell height = %f cm\n" % (n_azim, spacing, width, height)
            for i in range(0, m):
                print "Calculating parameters for track %d...\n" % (i)
                phi = (2 * math.pi / n_azim) * (i + 0.5)
                print "phi = %f" % (phi)
                nx = math.floor((width / spacing) * math.cos(phi)) + 1
                print "nx = %f" % (nx)
                ny = math.floor((height / spacing) * math.sin(phi)) + 1
                print "ny = %f " % (ny)
                phi_eff[i] = math.atan((height * nx)/(width * ny))
                print "phi_eff = %f " % (phi_eff[i])
                spacing_eff[i] = (width / nx) * math.sin(phi_eff[i])
                print "t_eff = %f cm\n" % (spacing_eff[i])
        else:
            print "Error: number of azimuthal angles desired must be a multiple of 4."
        return

    def azim_quadrature():


        """Function for calculating the weights for each azimuthal angle"""

        math = 1

        return