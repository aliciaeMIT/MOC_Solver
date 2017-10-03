import math


class TrackGeneration():
    def __init__(self, n_azim, spacing, width, height, m):
        self.n_azim = n_azim
        self.spacing = spacing
        self.width = width
        self.height = height
        self.m = m

    def phi(self, i):
        """computes azimuthal angle based on user desired inputs
        """
        return ((2 * math.pi / self.n_azim) * (i + 0.5))

    def nx(self, i):
        """ computes number of tracks for a given azimuthal angle, on y axis
        """
        return (math.floor((self.width / self.spacing) * math.cos(self.phi(i))) + 1)

    def ny(self, i):
        """computes number of tracks for given azimuthal angle on x axis
        """
        return (math.floor((self.height / self.spacing) * math.sin(self.phi(i))) + 1)

    def phi_eff(self, i):
        """computes effective azimuthal angle to guarantee cyclic tracks
        """
        return (math.atan((self.height * self.nx(i)) / (self.width * self.ny(i))))

    def spacing_eff(self, i):
        """ computation of effective ray spacing for each azimuthal angle
        """
        return ((self.width / self.nx(i)) * math.sin(self.phi_eff(i)))

    def azim_weight(self, i):
        """computation of azimuthal angle quadrature set, based on fraction of angular space of each angle.
        """
        if (i == 0):
            return ( (((self.phi_eff(i+1) - self.phi_eff(i)) / 2) + self.phi_eff(i)) / (2 * math.pi) )
        elif (i < (self.m - 1)):
            return ( (((self.phi_eff(i+1) - self.phi_eff(i)) / 2) + ((self.phi_eff(i) - self.phi_eff(i-1))/ 2)) / (2 * math.pi) )
        else:
            return ((2 * math.pi - self.phi_eff(i) + (self.phi_eff(i) - self.phi_eff(i-1)) / 2)/ (2 * math.pi))

    def polar_weight(self, n_p):
        """determines polar quadrature weights and angles using Tabuchi and Yamamoto (TY) set;
        number of polar divisions can be specified (2 or 3)
        returns ([polar weight per division], [sin-theta_polar per division])
        """
        if n_p == 2:
            wp = [0.212854 , 0.787146]
            sintheta_p = [0.363900 , 0.899900]
            return wp, sintheta_p
        elif n_p == 3:
            wp = [0.046233, 0.283619, 0.670148]
            sintheta_p = [0.166648, 0.537707, 0.932954]
            return wp, sintheta_p
        else:
            print "Error: must use 2 or 3 polar divisions for polar quadrature."
            return