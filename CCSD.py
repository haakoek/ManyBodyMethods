from numpy import *

def delta(i,j):
	if(i == j):
		return 1
	else:
		return 0
		
def QRPS(q,r,p,s,w):
	
	key1 = str(q+1)+str(r+1) + str(p+1) + str(s+1)
	key2 = str(q+1)+str(r+1) + str(s+1) + str(p+1)
	
	if((w.has_key(key1) == True) and (w.has_key(key2) == True)):
		qrps = w[key1]-0.5*w[key2]
	elif ((w.has_key(key1) == True) and (w.has_key(key2) == False)):
		qrps = w[key1]
	elif ((w.has_key(key1) == False) and (w.has_key(key2) == True)):
		qrps = -0.5*w[key2]
	else:
		qrps = 0.0
	return qrps

def energy(spatial_index):
	if spatial_index == 1:
		return 1.0
	elif ((spatial_index == 2) | (spatial_index == 3)):
		return 2.0
	elif ((spatial_index == 4) | (spatial_index == 5) | (spatial_index == 6)):
		return 3.0

def computeFockMatrix(N,L,w):
	F = zeros((L/2,L/2))
	for p in range(0,L/2):
		for q in range(0,L/2):
			F[p,q] = delta(p+1,q+1)*energy(q+1)
			for i in range(0,L/2):
				F[p,q] += QRPS(p+1,i+1,q+1,i+1,w)			
	return F	

inFile = open('coulomb.dat','r')
w = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key] = val

L = 12
N = 2
def initialize(N,L,w):
	F_init = computeFockMatrix(N,L,w)

	t1_old = zeros(N*(L/2-N))
	k = 0

	for i in range(0,N):
		for a in range(N,L/2):
			Dia  = F_init[i,i] - F_init[a,a]
			t_ia = F_init[a,i]/Dia
			t1_old[k] = t_ia
			k += 1

	k = 0
	t2_old = zeros(N**2 * (L/2-N)**2)		
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,L/2):
				for b in range(N,L/2):
					D_ij_ab = F_init[i,i] + F_init[j,j] - F_init[a,a] - F_init[b,b]
					t2_old[k] = QRPS(a+1,b+1,i+1,j+1,w)/D_ij_ab
					k += 1
	return F_init, t1_old, t2_old
				
F, t1_old, t2_old = initialize(N,L,w)
print t2_old
















