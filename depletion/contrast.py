import numpy as np
import matplotlib.pyplot as pl

# np.save("tb.npy",data)
def plot_proj(filename,second=2.0, color='k'):
	
	pot = np.loadtxt(filename)

	L = 4.0
	N = 10

	small = np.exp(-50)
	bins = np.linspace(0, L, N+1)

	# find the corresponding bin
	secondbin =np.digitize(second,bins)

	# pl.imshow(pot[here*N:here*(N+1),:])
	# pl.colorbar()
	# pl.show()
	for k in range(10):
		pl.plot(bins[:-1], pot[k*N:(k+1)*N,secondbin], label=filename, color=color,alpha=0.5)



plot_proj("Pe10pot.txt",color="b")
plot_proj("Pe60pot.txt",color="r")
# pl.legend()
pl.show()

