import numpy as np
from numpy import array
from pprint import pprint

class pt(object):
    def __init__(self, y=0, z=0):
        self.y = y
        self.z = z
    def show(self):
        print('pt = ', self.y, ',', self.z)
H = 8.
B = 8.
NFz = 8
NFy = 8
z_st = H/2*(1/NFz - 1)
y_st = B/2*(1/NFy - 1)
dz = H/NFz
dy = B/NFy
fiber_array = [[pt() for nfy in range(NFy)] for nfz in range(NFz)]

for nfz in range(NFz):
    for nfy in range(NFy):
        fiber_array_y = []
        y = y_st + B/NFy*nfy
        z = z_st + H/NFz*nfz
        fiber_array[nfy][nfz] = pt(y, z)
        
fiber_array[0][0].show()
