from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
import numpy as np
import pylab as pl
import sys
from pathlib import Path
import pickle

from collections import Counter
import sys
from stringato import extract_floats as ef

filename = sys.argv[1]

with open(filename, 'rb') as fin:
  N = int(fin.readline())
rho=float(ef(filename)[-1])
print("Rho", rho)
print("N", N)
L = (N/rho*0.5)**(1./3)

print("L",L)


node = import_file(filename,multiple_frames=True,columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z"])


def set_cell(frame, data):
    with data.cell_:
        data.cell_[:,0] = [L, 0., 0.]
        data.cell_[:,1] = [0., L, 0.]
        data.cell_[:,2] = [0., 0., 2*L]
        # #cell origin
        data.cell_[:,3] = [-L/2,  -L/2  ,  -L]
        # #set periodic boundary conditions
        data.cell_.pbc = (True, True, True)


node.modifiers.append(set_cell)

# modifier = WrapPeriodicImagesModifier()
# node.modifiers.append(modifier)

create_bonds_mod = VoronoiAnalysisModifier(compute_indices=False,generate_bonds=False)
node.modifiers.append(create_bonds_mod)


every = 1
hd = []
rhos= []
for frame in range(500, node.source.num_frames,every):
  data = node.compute(frame)
  vol = data.particles['Atomic Volume'].array
  rho = 1/vol
  rhos= np.concatenate((rhos,rho))
print (len(rhos)/N)
H,e=np.histogram(rhos, bins=35,density=True)
centre = e[:-1]+(e[1]-e[0])/2.0
np.savetxt(filename+".vol.hist.txt", list(zip(centre,H)))
# pl.show()
