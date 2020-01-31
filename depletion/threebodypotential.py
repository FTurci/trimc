import numpy as np
import h5py
import sys
import matplotlib.pyplot as pl
from stringato import extract_floats as ef
import numpy as np
from scipy.spatial.distance import pdist
from numba  import jit
import os

@jit(nopython=True)
def get_dists(box):
    points = np.random.uniform(0,1,size=(3,3))
    for k in range(3):
        points[:,k] *=box[k] 

    dists = np.zeros(3)
    count=0
    for i in range(2):
        for j in range(i+1,3):
            d=np.absolute(points[i]-points[j])
            d[d>L/2]-=L
            dist = np.linalg.norm(d)
            dists[count] = dist
            count+=1
    return sorted(dists)

@jit(nopython=True)
def repeat(nrepeat,box):
    data=np.zeros((nrepeat,3))
    for k in range(nrepeat):
        data[k]=get_dists(box)
    return data


def get_probability(values,bins):
    H,e = np.histogramdd(values, bins=(bins,bins,bins), normed=True)

    H[H==0]=H.min()
    return H

# filename = sys.argv[1]
# Pe = ef(filename)[1]
# L = float(ef(filename)[0])

Pe = 10.0
filename = f"store-L5.0Pe{Pe}.hdf5"

L = 5.0
N = 100

box = [L,L,L]

idealfile= "ideal-L%s.hdf5"%L
if os.path.isfile(idealfile):
    print("File already exists")
    with h5py.File(idealfile, 'r') as fin:
        ideal = fin['dataset'][()]
else:    
    ideal=repeat(100000000,box)#100000000
    print  ("Create ideal")
    with h5py.File(idealfile, 'w') as fout:
        fout.create_dataset('dataset',data=ideal)


# print (ideal)
small = np.exp(-50)
bins = np.linspace(0, L, N+1)
cutoff =2.5

ideal_pb = get_probability(ideal,bins)

from scipy.ndimage.filters import median_filter
with h5py.File(filename, 'r') as fin:
    # print(fin.keys())
    
    H = np.zeros((N,N,N))
    count= 0
    print(fin.keys())
    for key,value in fin.items():
        v =np.array(value)

        # if =="run2" or key=="run3" or key=="run4":
        if len(v)>1:
            print(v.shape)
            print("key",key)
            values = v

            # if (np.all((values[:,0]<=values[:,1])* (values[:,1]<=values[:,2]))):
            #     print("Order is ok")
            # else:
            #     print("Order is wrong")
            #     exit()


            HH,e = np.histogramdd(values, bins=(bins,bins,bins), normed=True)
            H += HH
            count+=1
    H/=count
    low_stat = H<H.max()/80
    print ("The ratio between the lowest and highest probabilities in the active sample is", H[H>0].min()/H.max())
    potential = -np.log(H/ideal_pb)
    # smoothen
    potential[np.logical_not(np.isfinite(potential))]=0
    
    # potential= median_filter(potential, size=4)
    
    potential[low_stat]=0

    with open(f"../Interactions/threebodypotential-L5-Pe{Pe}.txt", 'w') as fout:
        for s in range(potential.shape[2]):
            np.savetxt(fout,potential[:,:,s])
    np.save(f"../Interactions/threebodypotential-L5-Pe{Pe}.lowstat.npy",low_stat)
                
