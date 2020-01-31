import numpy as np
import h5py
import sys
import matplotlib.pyplot as pl
from stringato import extract_floats as ef
import numpy as np
from scipy.spatial.distance import pdist
from numba  import jit
import os
from scipy.interpolate import UnivariateSpline
from mayavi import mlab
import ipyvolume as ipv
from scipy.ndimage.filters import median_filter

def reassign(b,shape):
    dummy = np.zeros(shape)
    s0 = shape[0]
    for s in range(shape[2]):
        dummy[:,:,s] = b[s*s0:(s+1)*s0,:]
    return dummy

Pe =60.0
L = 5
N = 100
rmin=1.0

potential = np.loadtxt(f"../Interactions/threebodypotential-L5-Pe{Pe}.txt")

depletion = np.loadtxt(f"../Interactions/depletion-Pe{Pe}L4.0.txt")



tb = UnivariateSpline(depletion[:,0], depletion[:,1],s=0.01)
pl.plot(depletion[:,0], depletion[:,1])
pl.plot(depletion[:,0], tb(depletion[:,0]))
x = depletion[:,0]
pl.plot(x, 4*(1/x**12-1./x**6))
pl.show()
def twobody(r,rcut=2.5):
	return np.piecewise(r,[r<0.1, r>rcut, ],[0,0,tb])

def twobody_triangle(r1,r2,r3,rcut=2.5):
	U= twobody(r1)+twobody(r2)+twobody(r3)

	forbidden = (r1>r2)+(r3<r2)+(r3<r1)+(r1<rmin)+(r2<rmin)+(r3<rmin)#+(r1>rcut)+(r2>rcut)+(r3>rcut)
	U[forbidden]=0
	return U

def sticky_triangle(r1,r2,r3, epsilon,):
	U= (r1<2.5)*(r2<2.5)*(r3<2.5)*epsilon
	forbidden = (r1>r2)+(r3<r2)+(r3<r1)+(r1<rmin)+(r2<rmin)+(r3<rmin)+(r1>2.5)+(r2>2.5)+(r3>2.5)
	U[forbidden]=0
	return U

low_stat = np.load(f"../Interactions/threebodypotential-L5-Pe{Pe}.lowstat.npy",)
potential = reassign(potential,[N,N,N])
bins = np.linspace(0, L, N+1)
delta = bins[1]-bins[0]
cnt = bins[:-1]+(delta)/2


R1,R2,R3 = np.meshgrid(cnt,cnt,cnt, indexing="ij")

twopotential =  twobody_triangle(R1,R2,R3)

oldp = np.copy(potential)
zero = oldp==0

def condition(r1,r2,r3):
    forbidden = (r1>r2)+(r3<r2)+(r3<r1)+(r1<rmin)+(r2<rmin)+(r3<rmin)#+(r1>2.2)+(r2>2.2)+(r3>2.2)
    return forbidden 

forbidden = condition(R1,R2,R3)

potential = oldp-twopotential

potential= median_filter(potential, size=4)
# potential[forbidden]=0
potential[low_stat]=0
print(potential.max(), potential.min())
# potential[R2>5.0]=0
# potential[R3>5.0]=0
# potential[potential>0]=0
# potential = sticky_triangle(R1,R2,R3,-0.01)
with open(f"../Interactions/threebodypotential-minustwo-L5-Pe{Pe}.txt", 'w') as fout:
    for s in range(potential.shape[2]):
        np.savetxt(fout,potential[:,:,s])

fig, ax=pl.subplots(1,4, figsize=(8,3))
im=ax[0].imshow((potential).mean(axis=0).T, origin="left", extent=[cnt[0],cnt[-1],cnt[0],cnt[-1]])
# pl.xlabel("$r_1$"),pl.ylabel("$r_2$")
ax[0].set_xlabel("$r_2$"),ax[0].set_ylabel("$r_3$")
fig.colorbar(im, ax=ax[0])
im=ax[1].imshow((potential).mean(axis=1).T, origin="left", extent=[cnt[0],cnt[-1],cnt[0],cnt[-1]])
# pl.xlabel("$r_1$"),pl.ylabel("$r_2$")
ax[1].set_xlabel("$r_1$"),ax[0].set_ylabel("$r_3$")
fig.colorbar(im, ax=ax[1])

# pl.figure()
ax[3].hist(potential[oldp!=0].ravel(),bins=100)

# pl.imshow((potential).mean(axis=0).T, origin="left", extent=[cnt[0],cnt[-1],cnt[0],cnt[-1]])
# pl.xlabel("$r_2$"),pl.ylabel("$r_3$")
# pl.colorbar()

print("max  min", potential.max(),potential.min())


# energy of isosceles triangles at fixed angle 
n = 100
# angle = 65/180.*np.pi
# r1 = np.linspace(1,2.0,100)
# r2 = r1
# r3 = np.sqrt( r1**2+r2**2-2*r1*r2*np.cos(angle))

# angle = 45/180.*np.pi
# r2 = r3
# r1 = np.sqrt( r3**2+r2**2-2*r3*r2*np.cos(angle))

r1 = np.ones(n)*1.5
r2 = np.linspace(0.95,3.0,n)
r3 = np.ones(n)*2.

def f(a,b,c):
	i1 = int(a/delta)
	i2 = int(b/delta)
	i3 = int(c/delta)
	return potential[i1,i2,i3]


def ff(r1,r2,r3):
	u = np.zeros(len(r1))
	for i in range(len(r1)):
		a=r1[i]
		b=r2[i]
		c=r3[i]
		i1 = int(a/delta)
		i2 = int(b/delta)
		i3 = int(c/delta)
		u[i]=potential[i1,i2,i3]
		print(a,b,c,u[i], oldp[i1,i2,i3])
	return u
# ff =np.vectorize(f)

# pl.figure()
# pl.plot(r1,r3)

ax[2].plot(r2, twobody(r1)+twobody(r2)+twobody(r3)+ff(r1,r2,r3),'-o', label="all")
ax[2].plot(r2, twobody(r1)+twobody(r2)+twobody(r3), label="2")
ax[2].plot(r2, ff(r1,r2,r3), label="3")
ax[2].legend()
# for i in range(n):
	# print (r1[i],r2[i],r3[i], ff(r1[i],r2[i],r3[i]))
# pl.plot(r1, twobody(r1)+twobody(r2)+twobody(r3)+ff(r1,r2,r3),'-o', label="all")
# pl.plot(r1, twobody(r1)+twobody(r2)+twobody(r3), label="2")
# pl.plot(r1, ff(r1,r2,r3), label="3")
pl.show()

