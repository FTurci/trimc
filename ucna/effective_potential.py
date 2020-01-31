import numpy as np


def lj(r):
	return 4*(1/r**12-r**6)
def wca(r):
	u = np.zeros(len(r))
	u[r<2**(1./6.)]=lj(r)-lj(2**(1./6.))
	return u

d = 1.0 (scale)

delta  = 3 (dimensions)

# tau_a = read_input("peristence time")

# gamma = read_input("friction coeff")

# Da = 
tau_tilde= tau_a/gamma

