import numpy as np

def delta(i,j):
	if(i == j):
		return 1.0
	else:
		return 0

def SingleParticleEnergy(x):
	if (x == 0):
		return 1.0

	elif (x == 1 or x == 2):
		return 2.0

	elif (x == 3 or x == 4 or x == 5):
		return 3.0

	elif (x == 6 or x == 7 or x == 8 or x == 9):
		return 4.0

def map_index(x):
	
	if(x == 0 or x == 1):
		return 0

	elif(x == 2 or x == 3):
		return 1

	elif (x == 4 or x == 5):
		return 2

	elif (x == 6 or x == 7):
		return 3

	elif (x == 8 or x == 9):
		return 4

	elif (x == 10 or x == 11):
		return 5

	elif (x == 12 or x == 13):
		return 6

	elif (x == 14 or x == 15):
		return 7

	elif (x == 16 or x == 17):
		return 8

	elif (x == 18 or x == 19):
		return 9

def spin(x):
	if(   (x == 0) | (x == 2) | (x == 4) | (x == 6) | (x == 8) | (x == 10) ):
		return -1
	elif( (x == 1) | (x == 3) | (x == 5) | (x == 7) | (x == 9) | (x == 11) ):
		return 1


def QRPS(q,r,p,s,w):
	
	key1 = str(map_index(q))+str(map_index(r)) + str(map_index(p)) + str(map_index(s))
	key2 = str(map_index(q))+str(map_index(r)) + str(map_index(s)) + str(map_index(p))
	
	if((w.has_key(key1) == True) and (w.has_key(key2) == True)):
		qrps = w[key1]-0.5*w[key2]
	elif ((w.has_key(key1) == True) and (w.has_key(key2) == False)):
		qrps = w[key1]
	elif ((w.has_key(key1) == False) and (w.has_key(key2) == True)):
		qrps = -0.5*w[key2]
	else:
		qrps = 0.0
		
	return qrps

def QRPS2(q,r,p,s,w):
	
	key1 = str(map_index(q)) +str(map_index(r)) + str(map_index(p)) + str(map_index(s))
	key2 = str(map_index(q)) +str(map_index(r)) + str(map_index(s)) + str(map_index(p))

	if((w.has_key(key1) == True) and (w.has_key(key2) == True)):
		qrps = delta(spin(q),spin(p))*delta(spin(r),spin(s))*w[key1]-delta(spin(q),spin(s))*delta(spin(r),spin(p))*w[key2]
	elif ((w.has_key(key1) == True) and (w.has_key(key2) == False)):
		qrps = delta(spin(q),spin(p))*delta(spin(r),spin(s))*w[key1]
	elif ((w.has_key(key1) == False) and (w.has_key(key2) == True)):
		qrps = -delta(spin(q),spin(s))*delta(spin(r),spin(p))*w[key2]
	else:
		qrps = 0.0
		
	return qrps

def QRPS3(q,r,p,s,w):
	key = str(map_index(q)) +str(map_index(r)) + str(map_index(p)) + str(map_index(s))
	if(w.has_key(key) == True):
		return w[key]
	else:
		return 0

################################################################################
def computeEref(N,w):
	Eref = 0.0
	for i in range(0,N):
		Eref += SingleParticleEnergy(map_index(i))
		for j in range(0,N):
			Eref += 0.5*QRPS2(i,j,i,j,w) #( QRPS(map_index(i),map_index(j),map_index(i),map_index(j),w) )
	return Eref

def computeFockMatrix2(N,L,w):
	F = np.zeros( (L,L) )
	for p in range(0,L):
		F[p,p] = SingleParticleEnergy(map_index(p))
		for q in range(0,L):
			for i in range(0,N):
				F[p,q] += QRPS2(p,i,q,i,w)
				#F[p,q] += QRPS3(p,i,q,i,w)*delta(spin(p),spin(q))*delta(spin(i),spin(i))#QRPS2(p,i,q,i,w)
				#F[p,q] -= QRPS3(p,i,i,q,w)*delta(spin(p),spin(i))*delta(spin(i),spin(q))	
	return F

def initialize(N,L,w,F):
	
	#F_init = computeFockMatrix(N,L,w)
	tsingles = np.zeros((L-N,N))
	
	for a in range(N,L):
		for i in range(0,N):
			Dia  = F[i,i] - F[a,a]
			tia  = F[i,a]/Dia
			tsingles[a-N,i] = tia
	
	tdoubles = np.zeros( (L-N,L-N,N,N) )

	for a in range(N,L):
		for b in range(N,L):
			for i in range(0,N):
				for j in range(0,N):
					Dabij = F[i,i]+F[j,j] - F[a,a] - F[b,b]
					tdoubles[a-N,b-N,i,j] = QRPS2(a,b,i,j,w)/Dabij
	return tsingles, tdoubles

def ECCSD(E0,t1,t2,N,L,w,F):
	
	ECCSD = 0.0
	
	for i in range(0,N):
		for a in range(N,L):
			ECCSD += F[i,a]*t1[a-N,i]

	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,L):
				for b in range(N,L):
					tmp1 = QRPS2(i,j,a,b,w)
					tmp1 *= 0.25*t2[a-N,b-N,i,j]

					tmp2 = QRPS2(i,j,a,b,w)
					tmp2 *= 0.5*t1[a-N,i]*t1[b-N,j]

					ECCSD += tmp1 + tmp2

	return ECCSD

def computeT2Amplitudes(N,L,w,t1,t2,F):
	
	tdoubles = np.zeros( (L-N,L-N,N,N) )

	for i in range(0,N):
		for j in range(0,N):
			for a in range(N,L):
				for b in range(N,L):
					
					tdoubles[a-N,b-N,i,j] = QRPS2(i,j,a,b,w)

					for k in range(0,N):
						for l in range(0,N):
							tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(i,j,k,l,w)*t2[a-N,b-N,k,l]

					for c in range(N,L):
						for d in range(N,L):
							tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(c,d,a,b,w)*t2[c-N,d-N,i,j]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									tdoubles[a-N,b-N,i,j] += 0.25*QRPS2(c,d,k,l,w)*t2[a-N,b-N,k,l]

					for k in range(0,N):
						if(k != j):
							tdoubles[a-N,b-N,i,j] -= F[k,j]*t2[a-N,b-N,i,k]
						if(k != i):
							tdoubles[a-N,b-N,i,j] += F[k,i]*t2[a-N,b-N,j,k]

					for c in range(N,L):
						if(c != b):
							tdoubles[a-N,b-N,i,j] += F[b,c]*t2[a-N,c-N,i,j]
						if(c != a):
							tdoubles[a-N,b-N,i,j] -= F[a,c]*t2[b-N,c-N,i,j]

					#P(ab) -terms
					for k in range(0,N):
						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += F[k,c]*t1[a-N,k]*t2[b-N,c-N,i,j] - F[k,c]*t1[b-N,k]*t2[a-N,c-N,i,j]

					for k in range(0,N):
						tdoubles[a-N,b-N,i,j] -= QRPS2(i,j,k,b,w)*t1[a-N,k] - QRPS2(i,j,k,a,w)*t1[b-N,k]

					for k in range(0,N):
						for l in range(0,N):
							tdoubles[a-N,b-N,i,j] += 0.5*(QRPS2(i,j,k,l,w)*(t1[a-N,k]*t1[b-N,l] - t1[b-N,k]*t1[a-N,l]))

					for k in range(0,N):
						for c in range(N,L):
							for d in range(N,L):
								  tdoubles[a-N,b-N,i,j] += QRPS2(c,d,k,a,w)*t1[c-N,k]*t2[d-N,b-N,i,j] - 0.5*QRPS2(c,d,k,b,w)*t1[a-N,k]*t2[c-N,d-N,i,j]
								  tdoubles[a-N,b-N,i,j] -= QRPS2(c,d,k,b,w)*t1[c-N,k]*t2[d-N,a-N,i,j] - 0.5*QRPS2(c,d,k,a,w)*t1[b-N,k]*t2[c-N,d-N,i,j]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									tdoubles[a-N,b-N,i,j] += QRPS2(c-N,d-N,k,l,w)*(0.25*t1[a-N,k]*t1[b-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[a-N,l]*t2[d-N,b-N,i,j] - 0.5*t2[a-N,c-N,i,j]*t2[b-N,d-N,k,l])
									tdoubles[a-N,b-N,i,j] -= QRPS2(c-N,d-N,k,l,w)*(0.25*t1[b-N,k]*t1[a-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[b-N,l]*t2[d-N,a-N,i,j] - 0.5*t2[b-N,c-N,i,j]*t2[a-N,d-N,k,l])

					#P(ij)-terms
					for c in range(N,L):
						tdoubles[a-N,b-N,i,j] += QRPS2(c,j,a,b,w)*t1[c-N,i] - QRPS2(c,i,a,b,w)*t1[c-N,j]

					for c in range(N,L):
						for d in range(N,L):
							tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(c,d,a,b,w)*t1[c-N,i]*t1[d-N,j] - 0.5*QRPS2(c,d,a,b,w)*t1[c-N,j]*t1[d-N,i]

					for k in range(0,N):
						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += F[k,c]*t1[c-N,i]*t2[a-N,b-N,j,k] - F[k,c]*t1[c-N,j]*t2[a-N,b-N,i,k]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(c,j,k,l,w)*t1[c-N,i]*t2[a-N,b-N,k,l] - QRPS2(c,i,k,l,w)*t1[c-N,k]*t2[a-N,b-N,l,j]
								tdoubles[a-N,b-N,i,j] -= 0.5*QRPS2(c,i,k,l,w)*t1[c-N,j]*t2[a-N,b-N,k,l] - QRPS2(c,j,k,l,w)*t1[c-N,k]*t2[a-N,b-N,l,i]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								for d in range(N,L):							 						
									tdoubles[a-N,b-N,i,j] += QRPS2(c,d,k,l,w)*(0.25*t1[c-N,i]*t1[d-N,j]*t2[a-N,b-N,k,l] - t1[c-N,k]*t1[d-N,i]*t2[a-N,b-N,l,j] - 0.5*t2[a-N,b-N,i,k]*t2[c-N,d-N,j,l])
									tdoubles[a-N,b-N,i,j] -= QRPS2(c,d,k,l,w)*(0.25*t1[c-N,j]*t1[d-N,i]*t2[a-N,b-N,k,l] - t1[c-N,k]*t1[d-N,j]*t2[a-N,b-N,l,i] - 0.5*t2[a-N,b-N,j,k]*t2[c-N,d-N,i,l])

					#P(ij)P(ab)-terms

					for k in range(0,N):
						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += QRPS2(c,j,k,b,w)*t2[a-N,c-N,i,k] - QRPS2(i,c,k,b,w)*t1[a-N,k]*t1[c-N,j]
							tdoubles[a-N,b-N,i,j] -= QRPS2(c,j,k,a,w)*t2[b-N,c-N,i,k] - QRPS2(i,c,k,a,w)*t1[b-N,k]*t1[c-N,j]
							tdoubles[a-N,b-N,i,j] -= QRPS2(c,i,k,b,w)*t2[a-N,c-N,j,k] - QRPS2(j,c,k,b,w)*t1[a-N,k]*t1[c-N,i]
							tdoubles[a-N,b-N,i,j] += QRPS2(c,i,k,a,w)*t2[b-N,c-N,j,k] - QRPS2(j,c,k,a,w)*t1[b-N,k]*t1[c-N,i]

					for k in range(0,N):
						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += QRPS2(d,c,a,k,w)*t1[d-N,i]*t2[b-N,c-N,j,k] - QRPS2(c,d,k,b,w)*t1[c-N,i]*t1[a-N,k]*t1[d-N,j]
								tdoubles[a-N,b-N,i,j] -= QRPS2(d,c,b,k,w)*t1[d-N,i]*t2[a-N,c-N,j,k] - QRPS2(c,d,k,a,w)*t1[c-N,i]*t1[b-N,k]*t1[d-N,j]
								tdoubles[a-N,b-N,i,j] -= QRPS2(d,c,a,k,w)*t1[d-N,j]*t2[b-N,c-N,i,k] - QRPS2(c,d,k,b,w)*t1[c-N,j]*t1[a-N,k]*t1[d-N,i]	
								tdoubles[a-N,b-N,i,j] += QRPS2(d,c,b,k,w)*t1[d-N,j]*t2[a-N,c-N,i,k] - QRPS2(c,d,k,a,w)*t1[c-N,j]*t1[b-N,k]*t1[d-N,i]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(c,j,k,l,w)*t1[c-N,i]*t1[a-N,k]*t1[b-N,l] + QRPS2(i,c,k,l,w)*t1[a-N,l]*t2[b-N,c-N,j,k]
								tdoubles[a-N,b-N,i,j] -= 0.5*QRPS2(c,j,k,l,w)*t1[c-N,i]*t1[b-N,k]*t1[a-N,l] + QRPS2(i,c,k,l,w)*t1[b-N,l]*t2[a-N,c-N,j,k]
								tdoubles[a-N,b-N,i,j] -= 0.5*QRPS2(c,i,k,l,w)*t1[c-N,j]*t1[a-N,k]*t1[b-N,l] + QRPS2(j,c,k,l,w)*t1[a-N,l]*t2[b-N,c-N,i,k]
								tdoubles[a-N,b-N,i,j] += 0.5*QRPS2(c,i,k,l,w)*t1[c-N,j]*t1[b-N,k]*t1[a-N,l] + QRPS2(j,c,k,l,w)*t1[b-N,l]*t2[a-N,c-N,i,k]

					for k in range(0,N):
						for l in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									tdoubles[a-N,b-N,i,j] += QRPS2(c,d,k,l,w)*(t1[c-N,i]*t1[b-N,l]*t2[a-N,d-N,k,j] + 0.25*t1[c-N,i]*t1[a-N,k]*t1[d-N,j]*t1[b-N,l] + 0.5*t2[a-N,c-N,i,k]*t2[d-N,b-N,l,j])
									tdoubles[a-N,b-N,i,j] -= QRPS2(c,d,k,l,w)*(t1[c-N,i]*t1[a-N,l]*t2[b-N,d-N,k,j] + 0.25*t1[c-N,i]*t1[b-N,k]*t1[d-N,j]*t1[a-N,l] + 0.5*t2[b-N,c-N,i,k]*t2[d-N,a-N,l,j])
									tdoubles[a-N,b-N,i,j] -= QRPS2(c,d,k,l,w)*(t1[c-N,j]*t1[b-N,l]*t2[a-N,d-N,k,i] + 0.25*t1[c-N,j]*t1[a-N,k]*t1[d-N,i]*t1[b-N,l] + 0.5*t2[a-N,c-N,j,k]*t2[d-N,b-N,l,i])
									tdoubles[a-N,b-N,i,j] += QRPS2(c,d,k,l,w)*(t1[c-N,j]*t1[a-N,l]*t2[b-N,d-N,k,i] + 0.25*t1[c-N,j]*t1[b-N,k]*t1[d-N,i]*t1[a-N,l] + 0.5*t2[b-N,c-N,j,k]*t2[d-N,a-N,l,i])

					Dabij = F[i,i]+F[j,j] - F[a,a] - F[b,b]
					tdoubles[a-N,b-N,i,j] /= Dabij
	return tdoubles

#####################################################################################
inFile = open('coulomb2.dat','r')
w  = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key]  = val

N = 2; L = 6
prec = 1e-8

Eref = computeEref(N,w)
print "<phi0 | H | phi0 > = %g" % Eref

F = computeFockMatrix2(N,L,w)
t1_old, t2_old = initialize(N,L,w,F)

Eold = ECCSD(Eref,t1_old,t2_old,N,L,w,F)
print Eref+Eold
iters = 1

for k in range(1,20):

	t2_new = computeT2Amplitudes(N,L,w,t1_old,t2_old,F)
	Enew = ECCSD(Eref,t1_old,t2_new,N,L,w,F)
	
	iters += 1
	print Eref+Enew, abs(Enew-Eold), iters
	t2_old = t2_new

	if(abs(Enew-Eold) < prec):
		break

	Eold = Enew

#For N=2, L=6, ECCSD = 3.152329