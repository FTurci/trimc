import numpy as np
import matplotlib.pyplot as pl
from scipy.interpolate import UnivariateSpline

def reassign(b,shape):
    dummy = np.zeros(shape)
    s0 = shape[0]
    for s in range(shape[2]):
        dummy[:,:,s] = b[s*s0:(s+1)*s0,:]
    return dummy

depletion = np.loadtxt("../Interactions/depletion-Pe60.0L4.0.txt")

tb = UnivariateSpline(depletion[:,0], depletion[:,1],s=0)

def twobody(r):
	return np.piecewise(r,[r<1.0, r>2.5, ],[0,0,tb])
def twobody_triangle(r1,r2,r3):
	U= twobody(r1)+twobody(r2)+twobody(r3)

	forbidden = (r1>r2)+(r3<r2)+(r3<r1)+(r1<1)+(r2<1)+(r2<1)
	U[forbidden]=0
	return U

data = np.loadtxt("potential.txt")
u = reassign(data, (10,10,10))
L=5.0
N=10
bins = np.linspace(0, L, N+1)
R1,R2,R3 = np.meshgrid(bins[:-1],bins[:-1],bins[:-1])
cnt = bins[:-1]+(bins[1]-bins[0])/2.


R1,R2,R3 = np.meshgrid(cnt,cnt,cnt)

twopotential =  twobody_triangle(R2,R1,R3)

def triplets (R1,R2,R3):
	return 

def W (r1,r2,r3):
    atol = cnt[1]-cnt[0]
    i1,i2,i3= int(r1/atol),int(r2/atol),int(r3/atol)
    print(i1,i2,i3)
    return u[i1,i2,i3]



for s in range(u.shape[2]):
	print(bins[s])
	print(twopotential[:,:,s])
	print("-----------")
	pl.imshow(twopotential[:,:,s].T, origin="lower")
	pl.colorbar()
	pl.show()

print(W(1.6, 1.1,1.05))
