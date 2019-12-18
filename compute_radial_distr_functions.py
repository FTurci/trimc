from ovito.io import *
from ovito.data import *
from ovito.modifiers import *
import numpy as np
import pylab as pl
import sys



filename = sys.argv[1]
print(filename)
node = import_file(filename,multiple_frames=True,columns = ["Particle Type", "Position.X", "Position.Y", "Position.Z"])

bins = 300
rdf = np.zeros((bins,2), float)




print("Number of frames", node.source.num_frames)
collected=0
step=1#int(node.source.num_frames)
N =2048
rho=0.1
L = (N/rho)**(1./3.)
def set_cell(frame, data):
    with data.cell_:
        data.cell_[:,0] = [L, 0., 0.]
        data.cell_[:,1] = [0., L, 0.]
        data.cell_[:,2] = [0., 0., L]
        # #cell origin
        data.cell_[:,3] = [0,  0  ,  0]
        # #set periodic boundary conditions
        data.cell_.pbc = (True, True, True)

node.modifiers.append(set_cell)
modifier = CoordinationAnalysisModifier(cutoff = 6.0, number_of_bins = bins)
node.modifiers.append(modifier)


for frame in range(200,node.source.num_frames,step ):
  print(frame,)
  data = node.compute(frame)
  # cell=data.cell.matrix
  # L=cell[0,0]
  rdf += data.series['coordination-rdf'].as_table()

  collected+=1

print("Number of sampled frames", collected)
rdf/=collected

np.savetxt(filename+".rdf.txt",rdf)

pl.plot(rdf[:,0], rdf[:,1])

pl.show()
