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
				F[p,q] += QRPS3(p,i,q,i,w)*delta(spin(p),spin(q))*delta(spin(i),spin(i))#QRPS2(p,i,q,i,w)
				F[p,q] -= QRPS3(p,i,i,q,w)*delta(spin(p),spin(i))*delta(spin(i),spin(q))	
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
					tdoubles[a-N,b-N,i,j] = QRPS2(a,b,i,j,w)/Dabij#QRPS3(a,b,i,j,w)*delta(spin(a),spin(i))*delta(spin(b),spin(j))
					#tdoubles[a-N,b-N,i,j] -= QRPS3(a,b,j,i,w)*delta(spin(a),spin(j))*delta(spin(b),spin(i))
					#tdoubles[a-N,b-N,i,j] /= Dabij

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

#####################################################################################
inFile = open('coulomb2.dat','r')
w  = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key]  = val

N = 2; L = 6

Eref = computeEref(N,w)
print Eref

F = computeFockMatrix2(N,L,w)
t1_old, t2_old = initialize(N,L,w,F)

Enew = ECCSD(Eref,t1_old,t2_old,N,L,w,F)

print Eref+Enew
