import numpy as np
import matplotlib.pyplot as pl

x,y,d = np.loadtxt("../histogram.txt", unpack=True)

pl.plot(x,d)
pl.xlim(0,2)
pl.ylim(0,2)
#pl.plot(x,np.pi*4*x**2/3.5**3)
pl.show()
