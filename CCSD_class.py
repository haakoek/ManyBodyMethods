import numpy as np

class CCSD:
	def __init__(self, N,L,w,oneBodyElements, precision=1e-8):
		self.holeStates   		= N
		self.BasisFunctions     = L
		self.particleStates     = L-N
		self.twoBodyIntegrals   = w
		self.oneBodyIntegrals 	= oneBodyElements
		self.F    				= self.computeFockMatrix()
		self.t1_old 			= np.zeros((L-N,N))
		self.t2_old 			= np.zeros( (L-N,L-N,N,N) )
		self.precision   		= precision
	
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

	def initialize(self):
		
		N = self.holeStates
		L = self.BasisFunctions

		tsingles = np.zeros((L-N,N))
		
		for a in range(N,L):
			for i in range(0,N):
				Dia  = self.F[i,i] - self.F[a,a]
				tia  = self.F[i,a]/Dia
				self.t1_old[a-N,i] = tia
		
		for a in range(N,L):
			for b in range(N,L):
				for i in range(0,N):
					for j in range(0,N):
						Dabij = self.F[i,i]+self.F[j,j] - self.F[a,a] - self.F[b,b]
						self.t2_old[a-N,b-N,i,j] = self.QRPS2(a,b,i,j)/Dabij

	def ECCSD(self,t1,t2):
		
		N = self.holeStates
		L = self.BasisFunctions

		ECCSD = 0.0
		
		for i in range(0,N):
			for a in range(N,L):
				ECCSD += self.F[i,a]*t1[a-N,i]

		for i in range(0,N):
			for j in range(0,N):
				for a in range(N,L):
					for b in range(N,L):
						tmp1 = self.QRPS2(i,j,a,b)
						tmp1 *= 0.25*t2[a-N,b-N,i,j]

						tmp2 = self.QRPS2(i,j,a,b)
						tmp2 *= 0.5*t1[a-N,i]*t1[b-N,j]

						ECCSD += tmp1 + tmp2

		return ECCSD

	def solve(self):
		
		#Compute <phi0| H | phi0>
		Eref = self.computeEref()
		#Initialize amplitudes
		self.initialize()
		
		 
		Eold = self.ECCSD(self.t1_old,self.t2_old)
		print Eref+Eold
		
		#Compute 1 iteration
		t1_new = self.computeT1Amplitudes(self.t1_old,self.t2_old)
		t2_new = self.computeT2Amplitudes(self.t1_old,self.t2_old)
		Enew   = self.ECCSD(t1_new,t2_new)
		
		iters = 1
		print Eref+Enew, abs(Enew-Eold),iters
		self.t1_old = t1_new
		self.t2_old = t2_new
		
		

		while(abs(Enew-Eold) > self.precision):

			Eold   = Enew
			t1_new = self.computeT1Amplitudes(self.t1_old,self.t2_old)
			t2_new = self.computeT2Amplitudes(self.t1_old,self.t2_old)
			Enew   = self.ECCSD(t1_new,t2_new)
			
			iters += 1
			print Eref+Enew, abs(Enew-Eold), iters 
			self.t1_old = t1_new
			self.t2_old = t2_new
			

	def computeT1Amplitudes(self,t1,t2):
		
		N = self.holeStates
		L = self.BasisFunctions

		tsingles = np.zeros((L-N,N))

		for i in range(0,N):
			for a in range(N,L):
				
				tsingles[a-N,i] = self.F[i,a]

				for c in range(N,L):
					if (c != a):
						tsingles[a-N,i] += self.F[a,c]*t1[c-N,i]
				for k in range(0,N):
					if(i != k):
						tsingles[a-N,i] -= self.F[k,i]*t1[a-N,k]

				for k in range(0,N):
					for c in range(N,L):
						tsingles[a-N,i] += self.QRPS2(c,i,k,a)*t1[c-N,k]
						tsingles[a-N,i] += self.F[k,c]*t2[a-N,c-N,i,k]

				for k in range(0,N):
					for c in range(N,L):
						for d in range(N,L):
							tsingles[a-N,i] += 0.5*self.QRPS2(c,d,k,a)*t2[c-N,d-N,k,i]
							tsingles[a-N,i] += self.QRPS2(c,d,k,a)*t1[c-N,k]*t1[d-N,i]

				for k in range(0,N):
					for l in range(0,N):
						for c in range(N,L):
							tsingles[a-N,i] -= self.QRPS2(c,i,k,l)*t2[c-N,a-N,k,l]
							tsingles[a-N,i] -= self.QRPS2(c,i,k,l)*t1[c-N,k]*t1[a-N,l]
				
				for k in range(0,N):
					for c in range(N,L):
						tsingles[a-N,i] -= self.F[k,c]*t1[c-N,i]*t1[a-N,k]

				for k in range(0,N):
					for l in range(0,N):
						for c in range(N,L):
							for d in range(N,L):
								tsingles[a-N,i] -= self.QRPS2(c,d,k,l)*t1[c-N,k]*t1[d-N,i]*t1[a-N,l]
								tsingles[a-N,i] += self.QRPS2(c,d,k,l)*t1[c-N,k]*t2[d-N,a-N,l,i]
								tsingles[a-N,i] -= 0.5*self.QRPS2(c,d,k,l)*t2[c-N,d-N,k,i]*t1[a-N,l]
								tsingles[a-N,i] -= 0.5*self.QRPS2(c,d,k,l)*t2[c-N,a-N,k,l]*t1[d-N,i]

				Dia  = self.F[i,i] - self.F[a,a]
				tsingles[a-N,i] /= Dia

		return tsingles

	def computeT2Amplitudes(self,t1,t2):
	
		N = self.holeStates
		L = self.BasisFunctions

		tdoubles = np.zeros( (L-N,L-N,N,N) )

		for i in range(0,N):
			for j in range(0,N):
				for a in range(N,L):
					for b in range(N,L):
						
						tdoubles[a-N,b-N,i,j] = self.QRPS2(i,j,a,b)

						for k in range(0,N):
							for l in range(0,N):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(i,j,k,l)*t2[a-N,b-N,k,l]

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,d,a,b)*t2[c-N,d-N,i,j]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										tdoubles[a-N,b-N,i,j] += 0.25*self.QRPS2(c,d,k,l)*t2[a-N,b-N,k,l]

						for k in range(0,N):
							if(k != j):
								tdoubles[a-N,b-N,i,j] -= self.F[k,j]*t2[a-N,b-N,i,k]
							if(k != i):
								tdoubles[a-N,b-N,i,j] += self.F[k,i]*t2[a-N,b-N,j,k]

						for c in range(N,L):
							if(c != b):
								tdoubles[a-N,b-N,i,j] += self.F[b,c]*t2[a-N,c-N,i,j]
							if(c != a):
								tdoubles[a-N,b-N,i,j] -= self.F[a,c]*t2[b-N,c-N,i,j]

						#P(ab) -terms
						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += self.F[k,c]*t1[a-N,k]*t2[b-N,c-N,i,j] - self.F[k,c]*t1[b-N,k]*t2[a-N,c-N,i,j]

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] -= self.QRPS2(i,j,k,b)*t1[a-N,k] - self.QRPS2(i,j,k,a)*t1[b-N,k]

						for k in range(0,N):
							for l in range(0,N):
								tdoubles[a-N,b-N,i,j] += 0.5*(self.QRPS2(i,j,k,l)*(t1[a-N,k]*t1[b-N,l] - t1[b-N,k]*t1[a-N,l]))

						for k in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									  tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,a)*t1[c-N,k]*t2[d-N,b-N,i,j] - 0.5*self.QRPS2(c,d,k,b)*t1[a-N,k]*t2[c-N,d-N,i,j]
									  tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,b)*t1[c-N,k]*t2[d-N,a-N,i,j] - 0.5*self.QRPS2(c,d,k,a)*t1[b-N,k]*t2[c-N,d-N,i,j]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										tdoubles[a-N,b-N,i,j] += self.QRPS2(c-N,d-N,k,l)*(0.25*t1[a-N,k]*t1[b-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[a-N,l]*t2[d-N,b-N,i,j] - 0.5*t2[a-N,c-N,i,j]*t2[b-N,d-N,k,l])
										tdoubles[a-N,b-N,i,j] -= self.QRPS2(c-N,d-N,k,l)*(0.25*t1[b-N,k]*t1[a-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[b-N,l]*t2[d-N,a-N,i,j] - 0.5*t2[b-N,c-N,i,j]*t2[a-N,d-N,k,l])

						#P(ij)-terms
						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += self.QRPS2(c,j,a,b)*t1[c-N,i] - self.QRPS2(c,i,a,b)*t1[c-N,j]

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,d,a,b)*t1[c-N,i]*t1[d-N,j] - 0.5*self.QRPS2(c,d,a,b)*t1[c-N,j]*t1[d-N,i]

						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += self.F[k,c]*t1[c-N,i]*t2[a-N,b-N,j,k] - self.F[k,c]*t1[c-N,j]*t2[a-N,b-N,i,k]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,j,k,l)*t1[c-N,i]*t2[a-N,b-N,k,l] - self.QRPS2(c,i,k,l)*t1[c-N,k]*t2[a-N,b-N,l,j]
									tdoubles[a-N,b-N,i,j] -= 0.5*self.QRPS2(c,i,k,l)*t1[c-N,j]*t2[a-N,b-N,k,l] - self.QRPS2(c,j,k,l)*t1[c-N,k]*t2[a-N,b-N,l,i]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):							 						
										tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,l)*(0.25*t1[c-N,i]*t1[d-N,j]*t2[a-N,b-N,k,l] - t1[c-N,k]*t1[d-N,i]*t2[a-N,b-N,l,j] - 0.5*t2[a-N,b-N,i,k]*t2[c-N,d-N,j,l])
										tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,l)*(0.25*t1[c-N,j]*t1[d-N,i]*t2[a-N,b-N,k,l] - t1[c-N,k]*t1[d-N,j]*t2[a-N,b-N,l,i] - 0.5*t2[a-N,b-N,j,k]*t2[c-N,d-N,i,l])

						#P(ij)P(ab)-terms

						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += self.QRPS2(c,j,k,b)*t2[a-N,c-N,i,k] - self.QRPS2(i,c,k,b)*t1[a-N,k]*t1[c-N,j]
								tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,j,k,a)*t2[b-N,c-N,i,k] - self.QRPS2(i,c,k,a)*t1[b-N,k]*t1[c-N,j]
								tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,i,k,b)*t2[a-N,c-N,j,k] - self.QRPS2(j,c,k,b)*t1[a-N,k]*t1[c-N,i]
								tdoubles[a-N,b-N,i,j] += self.QRPS2(c,i,k,a)*t2[b-N,c-N,j,k] - self.QRPS2(j,c,k,a)*t1[b-N,k]*t1[c-N,i]

						for k in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									tdoubles[a-N,b-N,i,j] += self.QRPS2(d,c,a,k)*t1[d-N,i]*t2[b-N,c-N,j,k] - self.QRPS2(c,d,k,b)*t1[c-N,i]*t1[a-N,k]*t1[d-N,j]
									tdoubles[a-N,b-N,i,j] -= self.QRPS2(d,c,b,k)*t1[d-N,i]*t2[a-N,c-N,j,k] - self.QRPS2(c,d,k,a)*t1[c-N,i]*t1[b-N,k]*t1[d-N,j]
									tdoubles[a-N,b-N,i,j] -= self.QRPS2(d,c,a,k)*t1[d-N,j]*t2[b-N,c-N,i,k] - self.QRPS2(c,d,k,b)*t1[c-N,j]*t1[a-N,k]*t1[d-N,i]	
									tdoubles[a-N,b-N,i,j] += self.QRPS2(d,c,b,k)*t1[d-N,j]*t2[a-N,c-N,i,k] - self.QRPS2(c,d,k,a)*t1[c-N,j]*t1[b-N,k]*t1[d-N,i]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,j,k,l)*t1[c-N,i]*t1[a-N,k]*t1[b-N,l] + self.QRPS2(i,c,k,l)*t1[a-N,l]*t2[b-N,c-N,j,k]
									tdoubles[a-N,b-N,i,j] -= 0.5*self.QRPS2(c,j,k,l)*t1[c-N,i]*t1[b-N,k]*t1[a-N,l] + self.QRPS2(i,c,k,l)*t1[b-N,l]*t2[a-N,c-N,j,k]
									tdoubles[a-N,b-N,i,j] -= 0.5*self.QRPS2(c,i,k,l)*t1[c-N,j]*t1[a-N,k]*t1[b-N,l] + self.QRPS2(j,c,k,l)*t1[a-N,l]*t2[b-N,c-N,i,k]
									tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,i,k,l)*t1[c-N,j]*t1[b-N,k]*t1[a-N,l] + self.QRPS2(j,c,k,l)*t1[b-N,l]*t2[a-N,c-N,i,k]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,l)*(t1[c-N,i]*t1[b-N,l]*t2[a-N,d-N,k,j] + 0.25*t1[c-N,i]*t1[a-N,k]*t1[d-N,j]*t1[b-N,l] + 0.5*t2[a-N,c-N,i,k]*t2[d-N,b-N,l,j])
										tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,l)*(t1[c-N,i]*t1[a-N,l]*t2[b-N,d-N,k,j] + 0.25*t1[c-N,i]*t1[b-N,k]*t1[d-N,j]*t1[a-N,l] + 0.5*t2[b-N,c-N,i,k]*t2[d-N,a-N,l,j])
										tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,l)*(t1[c-N,j]*t1[b-N,l]*t2[a-N,d-N,k,i] + 0.25*t1[c-N,j]*t1[a-N,k]*t1[d-N,i]*t1[b-N,l] + 0.5*t2[a-N,c-N,j,k]*t2[d-N,b-N,l,i])
										tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,l)*(t1[c-N,j]*t1[a-N,l]*t2[b-N,d-N,k,i] + 0.25*t1[c-N,j]*t1[b-N,k]*t1[d-N,i]*t1[a-N,l] + 0.5*t2[b-N,c-N,j,k]*t2[d-N,a-N,l,i])

						Dabij = self.F[i,i] + self.F[j,j] - self.F[a,a] - self.F[b,b]
						tdoubles[a-N,b-N,i,j] /= Dabij
		return tdoubles

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
ccsd_test.solve()