import numpy as np
import sys
import os
from stringato import extract_floats as ef

path = sys.argv[1]
npart = sys.argv[2]

Pe = ef(path)[-2]

L=ef(path)[-1]



data = np.loadtxt(path)

data[:,0]= data[:,0]+0.5*(data[1,0]-data[0,0])
np.savetxt("tabular.txt",data)

os.system(f"./main {L} {npart} 1 {Pe}")

#command=f"mv ../histogram-type1-L{L}npart3.txt ../effective-Pe{Pe}-{L}-npart3.txt"
#print(command)
#os.system(command)
