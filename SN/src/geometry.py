from __future__ import division
import numpy as np


class Geometry(object):
    def __init__(self, pitch, mesh_spacing, fuel_width, fuel, moderator):

        self.width = pitch
        self.mesh = mesh_spacing
        self.fw = fuel_width
        self.fuel = fuel
        self.moderator = moderator

    def setMesh(self):
        self.n_cells = int(self.width / self.mesh) #number of cells in x, y
        n_cells = self.n_cells
        self.cells = np.zeros((n_cells, n_cells), dtype=object)

        n_mod = int((self.width / 2 - self.fw / 2) / self.mesh)
        n_fuel = int(self.fw / self.mesh)

        #setup 2d array of cells
        for j in range(n_cells):
            for i in range(n_cells):
                self.cells[i][j] = Cell()
                cell = self.cells[i][j]
                if (i > n_mod and j > n_mod) and \
                        ((i <= (n_mod + n_fuel) and j <= (n_mod + n_fuel))):
                    cell.region = 'fuel'
                    cell.getMaterial(self.fuel)
                    print "set cell %d, %d to %s " % (i, j, cell.region)
                else:
                    cell.getMaterial(self.moderator)
                    print "set cell %d, %d to %s " % (i, j, cell.region)


class Cell(object):
    def __init__(self):

        self.region = 'moderator'
        self.flux = 0
        self.angular = np.zeros(4)

    def getMaterial(self, material):
        self.material = material

class Material(object):
    def __init__(self, region, q, xs):

        self.name = region
        self.q = q
        self.xs = xs