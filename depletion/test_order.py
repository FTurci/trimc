import numpy as np

def reassign(b,a):
	dummy = np.zeros(a.shape)
	s0 = a.shape[0]
	for s in range(a.shape[2]):
		dummy[:,:,s] = b[s*s0:(s+1)*s0,:]
	return dummy

a = np.random.randint(10,size=(4,3,2))

with open("mat.txt", 'w') as fout:
	for s in range(a.shape[2]):
		np.savetxt(fout,a[:,:,s]) 

b = np.loadtxt("mat.txt")

b = reassign(b,a)
print (a)
print("==")
print (b)

print(a==b)