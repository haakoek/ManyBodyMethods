import numpy as np

class CCSD:
	def __init__(self, N,L,w,oneBodyElements):
		self.holeStates   		= N
		self.BasisFunctions     = L
		self.particleStates     = L-N
		self.twoBodyIntegrals   = w
		self.oneBodyIntegrals 	= oneBodyElements
	
	def delta(self,i,j):
		if(i == j):
			return 1.0
		else:
			return 0

	def map_index(self,x):
	
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

	def spin(self,x):
		if(   (x == 0) | (x == 2) | (x == 4) | (x == 6) | (x == 8) | (x == 10) ):
			return -1
		elif( (x == 1) | (x == 3) | (x == 5) | (x == 7) | (x == 9) | (x == 11) ):
			return 1

	def QRPS2(self,q,r,p,s):
	
		key1 = str(self.map_index(q)) +str(self.map_index(r)) + str(self.map_index(p)) + str(self.map_index(s))
		key2 = str(self.map_index(q)) +str(self.map_index(r)) + str(self.map_index(s)) + str(self.map_index(p))

		if((self.twoBodyIntegrals.has_key(key1) == True) and (self.twoBodyIntegrals.has_key(key2) == True)):
			qrps = self.delta(self.spin(q),self.spin(p))*self.delta(self.spin(r),self.spin(s))*self.twoBodyIntegrals[key1]-self.delta(self.spin(q),self.spin(s))*self.delta(self.spin(r),self.spin(p))*self.twoBodyIntegrals[key2]
		elif ((self.twoBodyIntegrals.has_key(key1) == True) and (self.twoBodyIntegrals.has_key(key2) == False)):
			qrps = self.delta(self.spin(q),self.spin(p))*self.delta(self.spin(r),self.spin(s))*self.twoBodyIntegrals[key1]
		elif ((self.twoBodyIntegrals.has_key(key1) == False) and (self.twoBodyIntegrals.has_key(key2) == True)):
			qrps = -self.delta(self.spin(q),self.spin(s))*self.delta(self.spin(r),self.spin(p))*self.twoBodyIntegrals[key2]
		else:
			qrps = 0.0
			
		return qrps
	def computeEref(self):
		Eref = 0.0
		for i in range(0,N):
			Eref += self.oneBodyIntegrals[i]
			for j in range(0,N):
					Eref += 0.5*self.QRPS2(i,j,i,j) #( QRPS(map_index(i),map_index(j),map_index(i),map_index(j),w) )
		return Eref

	def computeFockMatrix(self):
		
		N = self.holeStates
		L = self.BasisFunctions

		F = np.zeros( (L,L) )
		for p in range(0,L):
			F[p,p] = self.oneBodyIntegrals[p]
			for q in range(0,L):
				for i in range(0,N):
					F[p,q] += self.QRPS2(p,i,q,i)	
		return F

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

def SingleParticleEnergy(x):
	if (x == 0):
		return 1.0

	elif (x == 1 or x == 2):
		return 2.0

	elif (x == 3 or x == 4 or x == 5):
		return 3.0

	elif (x == 6 or x == 7 or x == 8 or x == 9):
		return 4.0

#################################################################################################################################################
import sys
inFile = open('coulomb2.dat','r')
w  = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key]  = val

N = int(sys.argv[1]); L = int(sys.argv[2])
prec = 1e-8

oneBodyElements = np.zeros(L)

for i in range(0,L):
	oneBodyElements[i] = SingleParticleEnergy(map_index(i))

ccsd_test = CCSD(N,L,w,oneBodyElements)

Eref = ccsd_test.computeEref()
F    = ccsd_test.computeFockMatrix()

print F
print Eref