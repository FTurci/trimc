import matplotlib.pyplot as pl
import numpy  as np
import glob
from natsort import natsorted
from pronto import Colors
from stringato import extract_floats as ef


files = natsorted(glob.glob("*.vol.hist.txt"))
colors = Colors(len(files))
pl.figure(figsize=(1.6,1.2))
for i,f in enumerate(files):
	c,H = np.loadtxt(f, unpack=True)
	rho = ef(f)[-1]
	lb = rho
	if i ==0 :
		lb = r"$\rho = $"+rho
	pl.plot(c,H, label=lb, color= colors.palette[i])
pl.legend(loc='upper left', ncol=2, handlelength=0.3,fontsize=8,labelspacing=0.2, handletextpad=0.3)
pl.xlabel(r"$\rho_{loc}$")
pl.ylabel(r"$p(\rho_{loc})$")
pl.ylim(0,8)
pl.yticks([0,4,8.])
pl.savefig("figtwoBodyOnlyVoroRhos.pdf", transparent=True)
pl.show()