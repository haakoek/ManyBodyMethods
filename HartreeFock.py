from numpy import *
import copy

def energy(spatial_index):
	if spatial_index == 1:
		return 1.0
	elif ((spatial_index == 2) | (spatial_index == 3)):
		return 2.0
	elif ((spatial_index == 4) | (spatial_index == 5) | (spatial_index == 6)):
		return 3.0
	elif ( (spatial_index == 7) | (spatial_index == 8) | (spatial_index == 9) | (spatial_index == 10) ):
		return 4.0
		
def delta(i,j):
	if(i == j):
		return 1
	else:
		return 0

def density_matrix(U,L,N):
	D = zeros((L/2,L/2))
	for r in range(0,L/2):
		for s in range(0,L/2):
			D[s,r] = 0.0
			for j in range(0,N/2):
				D[s,r] += 2.0*U[s,j]*conjugate(U[r,j])		
	return D
	
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

def compute_FockMatrix(D,L,w):
	FD = zeros((L/2,L/2))
	for q in range(0,L/2):
		for p in range(0,L/2):
			FD[q,p] = delta(q+1,p+1)*energy(p+1)
			for r in range(0,L/2):
				for s in range(0,L/2):
					FD[q,p] += D[s,r]*QRPS(q,r,p,s,w)
	return FD
	

def diag_FD(FD):
	eps, U = linalg.eig(FD)
	eps, U = bubbleSort(eps,U)
	return eps, U
	
def bubbleSort(alist,eig_vec):
	
	sorted_alist = copy.copy(alist)
	sorted_eigvec = copy.copy(eig_vec)
	
	for passnum in range(len(sorted_alist)-1,0,-1):
		for i in range(passnum):
			if sorted_alist[i]>sorted_alist[i+1]:
				temp = sorted_alist[i]
				temp2 = copy.copy(sorted_eigvec[:,i])
				sorted_alist[i] = sorted_alist[i+1]
				sorted_alist[i+1] = temp
				sorted_eigvec[:,i] = copy.copy(sorted_eigvec[:,i+1])
				sorted_eigvec[:,i+1] = temp2			
	return sorted_alist, sorted_eigvec
	
def check_F_Hermitian(F,L):
	isHermitian = True
	for p in range(0,L/2):
		for q in range(0,L/2):
			if (F[p,q] != conjugate(F[q,p])):
				isHermitian = False
				break
	return isHermitian
	
def HF_Energy(eps, U, D, w, N, L):
	
	ERHF = 0.0
	for i in range(0,N/2):
		ERHF += eps[i]
	ERHF = 2.0*ERHF
	
	for i in range(0,N/2):
		for p in range(0,L/2):
			for q in range(0,L/2):
				tmp = 0.0
				for s in range(0,L/2):
					for r in range(0,L/2):
						tmp += D[s,r]*QRPS(q,r,p,s,w)
				ERHF -= conjugate(U[q,i])*tmp*U[p,i]
	return ERHF

def computeHartreeFockSolution(L,N,w):

	epsilon = 1*10**(-16)
	U = identity(L/2)
	D = density_matrix(U,L,N)
	FD = compute_FockMatrix(D,L,w)
	eps_old, U = diag_FD(FD)
	max_iters = 30
	ERHF = HF_Energy(eps_old,U,D,w,N,L)
	print "ERHF = %.16f, from initial guess." % ERHF

	for k in range(1,max_iters):
		D = density_matrix(U,L,N)
		FD = compute_FockMatrix(D,L,w)
		eps_new, U = diag_FD(FD)
		ERHF = HF_Energy(eps_new,U,D,w,N,L)
		#print eps_new
		print "ERHF = %.18f, |eps_new - eps_old| = %g" % (ERHF,max(abs(eps_new-eps_old)))	
		if(max(abs(eps_new-eps_old)) < epsilon):
			eps_old = eps_new		
			break	
		eps_old = eps_new

	return FD, U, eps_old, ERHF
