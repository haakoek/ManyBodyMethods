from numpy import *

def delta(i,j):
	if(i == j):
		return 1
	else:
		return 0
		
def QRPS(q,r,p,s,w):
	
	key1 = str(q)+str(r) + str(p) + str(s)
	key2 = str(q)+str(r) + str(s) + str(p)
	
	if((w.has_key(key1) == True) and (w.has_key(key2) == True)):
		qrps = w[key1]-0.5*w[key2]
	elif ((w.has_key(key1) == True) and (w.has_key(key2) == False)):
		qrps = w[key1]
	elif ((w.has_key(key1) == False) and (w.has_key(key2) == True)):
		qrps = -0.5*w[key2]
	else:
		qrps = 0.0
		
	return qrps

def QR_v_PS(q,r,p,s,w):
	key = str(q)+str(r)+str(p)+str(s)
	if(w.has_key(key)):
		return w[key]
	else:
		return 0.0

def energy(spatial_index):
	if spatial_index == 1:
		return 1.0
	elif ((spatial_index == 2) | (spatial_index == 3)):
		return 2.0
	elif ((spatial_index == 4) | (spatial_index == 5) | (spatial_index == 6)):
		return 3.0
	elif ( (spatial_index == 7) | (spatial_index == 8) | (spatial_index == 9) | (spatial_index == 10) ):
		return 4.0

def computeFockMatrix(N,L,w):
	F = zeros((L,L))
	for p in range(0,L):
		for q in range(0,L):
			F[p,q] = delta(p+1,q+1)*energy(q+1)
			for i in range(0,L):
				F[p,q] += QRPS(p,i,q,i,w)			
	return F

def initialize(N,L,w,F_init):
	
	#F_init = computeFockMatrix(N,L,w)

	t1_old = zeros((N, (L-N)))
	
	for i in range(0,N):
		for a in range(N,L):
			Dia  = F_init[i,i] - F_init[a,a]
			t_ia = F_init[a,i]/Dia
			t1_old[i,a-N] = t_ia
	
	k = 0
	t2_old = zeros((N,N,(L-N),(L-N)))		
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,L):
				for b in range(N,L):
					D_ij_ab = F_init[i,i] + F_init[j,j] - F_init[a,a] - F_init[b,b]
					t2_old[i,j,a-N,b-N] = QRPS(a,b,i,j,w)/D_ij_ab
	return t1_old, t2_old	

#Extremely naive implementation of coupled-cluster amplitude equations.

def ECCSD(E0,t1,t2,N,L,w,F):
	E_CCSD = E0
	
	for i in range(0,N):
		for a in range(N,L):
			E_CCSD += F[a,i]*t1[i,a-N]

	for a in range(N,L):
		for b in range(N,L):
			for i in range(0,N):
				for j in range(0,N):
					E_CCSD += 0.25*I(i,j,a,b,w)*t2[i,j,a-N,b-N]
					E_CCSD += 0.5*I(i,j,a,b,w)*t1[i,a-N]*t1[j,b-N]

	return E_CCSD

def I(i,j,a,b,w):
	return QRPS(a,b,i,j,w)

def computeT1Amplitudes(N,L,w,F,t1_old,t2_old):
	
	t1_new = zeros((N, (L-N)))
	
	for i in range(0,N):
		for a in range(N,L):
			
			t1_new[i,a-N] = F[a,i] # *D_ia
			
			for c in range(N,L):
				if(c != a):
					t1_new[i,a-N] += F[a,c]*t1_old[i,c-N]
			
			for k in range(0,N):
				if(i != k):
					t1_new[i,a-N] -= F[k,i]*t1_old[k,a-N]
			
			for k in range(0,N):
				for c in range(N,L):
					t1_new[i,a-N] += QRPS(c,i,k,a,w)*t1_old[k,c-N]
					t1_new[i,a-N] += F[k,c]*t2_old[i,k,a-N,c-N]
					t1_new[i,a-N] -= F[k,c]*t1_old[i,c-N]*t1_old[k,a-N]

			for k in range(0,N):
				for c in range(N,L):
					for d in range(N,L):
						t1_new[i,a-N] += 0.5*QRPS(c,d,k,a,w)*t2_old[k,i,c-N,d-N]
						t1_new[i,a-N] += QRPS(c,d,k,a,w)*t1_old[k,c-N]*t1_old[i,d-N]

			for k in range(0,N):
				for l in range(0,N):
					for c in range(N,L):
						t1_new[i,a-N] -= 0.5*QRPS(c,i,k,l,w)*t2_old[k,l,c-N,a-N]
						t1_new[i,a-N] -= QRPS(c,i,k,l,w)*t1_old[k,c-N]*t1_old[l,a-N]
			
			for k in range(0,N):
				for l in range(0,N):
					for c in range(N,L):
						for d in range(N,L):							
							t1_new[i,a-N] -= QRPS(c,d,k,l,w)*t1_old[k,c-N]*t1_old[i,d-N]*t1_old[l,a-N] 
							t1_new[i,a-N] += QRPS(c,d,k,l,w)*t1_old[k,c-N]*t2_old[l,i,d-N,a-N]
							t1_new[i,a-N] -= 0.5*QRPS(c,d,k,l,w)*t2_old[k,i,c-N,d-N]*t1_old[l,a-N]
							t1_new[i,a-N] -= 0.5*QRPS(c,d,k,l,w)*t2_old[k,l,c-N,a-N]*t1_old[i,d-N]
			
			D_ia = F[i,i]-F[a,a]			
			t1_new[i,a-N] = t1_new[i,a-N]/D_ia

	return t1_new

def computeT2Amplitudes(N,L,w,F,t1_old,t2_old):
	t2_new = zeros((N,N,(L-N),(L-N)))
	
	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,L):
				for b in range(N,L):
					
					D_ij_ab = F[i,i] + F[j,j] - F[a,a] - F[b,b]					
					t2_new[i,j,a-N,b-N] = I(a,b,i,j,w)

					for k in range(0,N):
						for l in range(0,N):
							t2_new[i,j,a-N,b-N] += 0.5*I(k,l,i,j,w)*t2_old[k,l,a-N,b-N]
							t2_new[i,j,a-N,b-N] += 0.5*(I(k,l,i,j,w)*t1_old[k,a-N]*t1_old[l,b-N] - I(k,l,i,j,w)*t1_old[k,b-N]*t1_old[l,a-N])

					for c in range(N,L):
						for d in range(N,L):
							t2_new[i,j,a-N,b-N] += 0.5*I(a,b,c,d,w)*t2_old[i,j,c-N,d-N]

					for k in range(0,N):
						for c in range(N,L):
							t2_new[i,j,a-N,b-N] += F[k,c]*t1_old[k,a-N]*t2_old[i,j,b-N,c-N] - F[k,c]*t1_old[k,b-N]*t2_old[i,j,a-N,c-N]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									t2_new[i,j,a-N,b-N] += 0.25*I(k,l,c,d,w)*t2_old[i,j,c-N,d-N]*t2_old[k,l,a-N,b-N]

					for k in range(0,N):
						if(k != j):						
							t2_new[i,j,a-N,b-N] -= F[k,j]*t2_old[i,k,a-N,b-N]
						if(k != i):
							t2_new[i,j,a-N,b-N] += F[k,i]*t2_old[j,k,a-N,b-N]

						t2_new[i,j,a-N,b-N] -= I(k,b,i,j,w)*t1_old[k,a-N] - I(k,a,i,j,w)*t1_old[k,b-N]

					for c in range(N,L):
						if(c != b):
							t2_new[i,j,a-N,b-N] += F[b,c]*t2_old[i,j,a-N,c-N]
						if(c != a):
							t2_new[i,j,a-N,b-N] -= F[a,c]*t2_old[i,j,b-N,c-N] 
					
					t2_new[i,j,a-N,b-N] = t2_new[i,j,a-N,b-N]/D_ij_ab  
	return t2_new

def computeEref():
	Eref = 0.0
	for i in range(0,N/2):
		Eref += 2*energy(i+1)
		for j in range(0,N/2):
			Eref += 0.5*( QRPS(i,j,i,j,w) )
	return Eref

from HartreeFock import *

#inFile = open('coulomb2.dat','r')
inFile = open('coulomb2.dat','r')
w  = {}
w2 = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key]  = 4*val
	w2[key] = val 

L = 12
N = 6

F, U, eps_old, ERHF = computeHartreeFockSolution(L,N,w2)

Eref = computeEref()
print Eref
F_init = computeFockMatrix(N/2,L/2,w2)
t1_old, t2_old = initialize(N/2,L/2,w2,F_init)
ECCSD = ECCSD(Eref,t1_old,t2_old,N/2,L/2,w2,F_init)
print ECCSD

"""
t1_old, t2_old = initialize(N/2,L,w,F)
Energy = ECCSD(ERHF,t1_old,t2_old,N/2,L,w,F)

print Energy
"""

"""
for i in range(1,20):
	t1_new = computeT1Amplitudes(N/2,L,w,F,t1_old,t2_old)
	t2_new = computeT2Amplitudes(N/2,L,w,F,t1_old,t2_old)
	Energy = ECCSD(ERHF,t1_new,t2_new,N/2,L,w,F)
	t1_old = t1_new
	t2_old = t2_new

print Energy
""" 

"""
t1_new = computeT1Amplitudes(N,L,w,F,t1_old,t2_old)
t2_new = computeT2Amplitudes(N,L,w,F,t1_old,t2_old)
 
Enew = ECCSD(0,t1_new,t2_new,N,L,w,F)
print Enew
"""


















