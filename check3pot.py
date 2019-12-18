import numpy as np
import matplotlib.pyplot as pl

data = np.loadtxt("tb.txt", usecols=[1,2,3,4])

# np.save("tb.npy",data)

pot = np.loadtxt("threebodypotential.txt")

L = 4.0
N = 10

small = np.exp(-50)
bins = np.linspace(0, L, N+1)

# fix third length to the average third length

third = data[:,2].mean()
second = data[:,1].mean()
# find the corresponding bin
thirdbin =np.digitize(third,bins)
secondbin =np.digitize(second,bins)

print(third, second)
# pl.imshow(pot[here*N:here*(N+1),:])
# pl.colorbar()
# pl.show()

pl.plot(bins[:-1]	, pot[thirdbin*N:(thirdbin+1)*N,secondbin])

tol = 0.1
selection = np.isclose(data[:,1],second,rtol=tol)*np.isclose(data[:,2],third,rtol=tol)
pl.plot(data[selection,0], data[selection,-1],'o', alpha=0.2)
pl.show()

