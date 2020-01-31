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

# filename = sys.argv[1]
N=864
rho=1.2
L = (N/rho)**(1./3)

print("L",L)


fig,ax =pl.subplots(2)

def high_rho(filename,Pe,nbins=21):
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

  
  every = 100
  hd = []
  for frame in range(0, node.source.num_frames,every):
    data = node.compute(frame)
    vol = data.particles['Atomic Volume'].array
    rho = 1/vol
    H,e = np.histogram(rho, bins=nbins,density=True)
    cntr = e[:-1]+(e[1]-e[0])/2
    # hd.append(cntr[cntr>0.8][H[cntr>0.8].argmax()])
    hd.append(H[cntr>0.8].sum()*(cntr[1]-cntr[0]))
  ax[0].plot(hd, label=str(Pe))
  
  ax[1].plot(cntr,H)
  np.savetxt(f"forDatagraph/Pe{Pe}.txt", list(zip(cntr,H)))
  # kde = 
  # rho_range = np.linspace(0.2,1.4,1000)
  # log_dens = kde.score_samples(rho_range)
  # pl.plot(rho_range, np.exp(log_dens))
  # ax[1].hist(rho, histtype="step", density=True, bins=nbins)

filenames={1.0:'rho1.2-Pe1.0.xyz',10.0:'rho1.2-Pe10.0.xyz',20.0:'rho1.2-Pe20.0.xyz',30.0:'rho1.2-Pe30.0.xyz',40.0:'rho1.2-Pe40.0.xyz',60.0:'rho1.2.xyz'}


for key,filename in filenames.items():
  high_rho(filename,key)
ax[0].legend(frameon=False, title="Pe", ncol=2)
ax[0].set_ylabel(r"$\int_{0.8}^{+\infty} p(\rho)d\rho$")
ax[0].set_xlabel("MC sweep/1000")
ax[1].set_xlabel(r"$\rho_{local}$")
pl.tight_layout()
pl.savefig("volumes.pdf")
pl.show()
