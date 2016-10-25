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

		self.W1 				= np.zeros((N,N,N,N))
		self.W2 				= np.zeros((L-N,N,N,N))
		self.W3 				= np.zeros((N,N,L-N,N))
		self.W4 				= np.zeros((L-N,N,N,L-N))
		self.F1 				= np.zeros((L-N,N))
		self.F2 				= np.zeros((N,N))
		self.F3 				= np.zeros((L-N,L-N))

	def updateIntermediates(self):

		for k in range(0,N):
			for l in range(0,N):
				for i in range(0,N):
					for j in range(0,N):
						self.W1[k,l,i,j] = self.computeW1(k,l,i,j)

		for a in range(N,L):
			for k in range(0,N):
				for i in range(0,N):
					for j in range(0,N):
						self.W2[a-N,k,i,j] = self.computeW2(a,k,i,j)
		
		for k in range(0,N):
			for l in range(0,N):
				for c in range(N,L):
					for i in range(0,N):
						self.W3[k,l,c-N,i] = self.computeW3(k,l,c,i)
		
		for c in range(N,L):
			for k in range(0,N):
				self.F1[c-N,k] = self.computeF1(c,k)

		for k in range(0,N):
			for i in range(0,N):
				self.F2[k,i] = self.computeF2(k,i)

		for a in range(N,L):
			for k in range(0,N):
				for i in range(0,N):
					for c in range(N,L):
						self.W4[a-N,k,i,c-N] = self.computeW4(a,k,i,c)

		for a in range(N,L):
			for c in range(N,L):
				self.F3[a-N,c-N] = self.computeF3(a,c)
		
	def updateIntermediatesCCD(self):

		for k in range(0,N):
			for l in range(0,N):
				for i in range(0,N):
					for j in range(0,N):
						self.W1[k,l,i,j] = self.computeW1CCD(k,l,i,j)

		for a in range(N,L):
			for k in range(0,N):
				for i in range(0,N):
					for j in range(0,N):
						self.W2[a-N,k,i,j] = self.computeW2CCD(a,k,i,j)
		
		for k in range(0,N):
			for l in range(0,N):
				for c in range(N,L):
					for i in range(0,N):
						self.W3[k,l,c-N,i] = self.computeW3CCD(k,l,c,i)
		
		for c in range(N,L):
			for k in range(0,N):
				self.F1[c-N,k] = self.computeF1CCD(c,k)

		for k in range(0,N):
			for i in range(0,N):
				self.F2[k,i] = self.computeF2CCD(k,i)

		for a in range(N,L):
			for k in range(0,N):
				for i in range(0,N):
					for c in range(N,L):
						self.W4[a-N,k,i,c-N] = self.computeW4CCD(a,k,i,c)


		for a in range(N,L):
			for c in range(N,L):
				self.F3[a-N,c-N] = self.computeF3CCD(a,c)

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
		if(   (x == 0) | (x == 2) | (x == 4) | (x == 6) | (x == 8) | (x == 10) | (x == 12) | (x == 14) | (x == 16) | (x == 18) ):
			return -1
		elif( (x == 1) | (x == 3) | (x == 5) | (x == 7) | (x == 9) | (x == 11) | (x == 13) | (x == 15) | (x == 17) | (x == 19) ):
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

	def ECCD(self,t2):

		N = self.holeStates
		L = self.BasisFunctions

		ECCD = 0.0
		
		for i in range(0,N):
			for j in range(0,N):
				for a in range(N,L):
					for b in range(N,L):
						tmp1 = self.QRPS2(i,j,a,b)
						tmp1 *= 0.25*t2[a-N,b-N,i,j]
						ECCD += tmp1

		return ECCD

	def solve(self):
		
		#Compute <phi0| H | phi0>
		Eref = self.computeEref()
		print Eref
		#Initialize amplitudes
		self.initialize()
		
		Eold = self.ECCSD(self.t1_old,self.t2_old)
		print Eref+Eold
		
		#Compute 1 iteration
		t1_new = self.computeT1Amplitudes(self.t1_old,self.t2_old)
		t2_new = self.computeT2Amplitudes(self.t1_old,self.t2_old)
		print t1_new

		Enew   = self.ECCSD(t1_new,t2_new)
		
		iters = 1
		print Eref+Enew, abs(Enew-Eold),iters
		self.t1_old = t1_new
		self.t2_old = t2_new
		
		"""	
		while(abs(Enew-Eold) > self.precision):

			Eold   = Enew
			t1_new = self.computeT1Amplitudes(self.t1_old,self.t2_old)
			t2_new = self.computeT2Amplitudes(self.t1_old,self.t2_old)
			Enew   = self.ECCSD(t1_new,t2_new)
			
			iters += 1
			print Eref+Enew, abs(Enew-Eold), iters 
			self.t1_old = t1_new
			self.t2_old = t2_new
		"""	
	
	def solveWithCCD(self):

		#Compute <phi0| H | phi0>
		Eref = self.computeEref()
		#Initialize amplitudes
		self.initialize()
		
		self.updateIntermediatesCCD()
		 
		Eold = self.ECCD(self.t2_old)
		print Eref+Eold
		
		#Compute 1 iteration
		
		t2_new = self.computeT2AmplitudesCCD(self.t2_old)
		Enew   = self.ECCD(t2_new)
		
		iters = 1
		print Eref+Enew, abs(Enew-Eold),iters
		self.t2_old = t2_new
				
		while(abs(Enew-Eold) > self.precision):

			Eold   = Enew
			t2_new = self.computeT2AmplitudesCCD(self.t2_old)
			Enew   = self.ECCD(t2_new)
			
			iters += 1
			print Eref+Enew, abs(Enew-Eold), iters 
			self.t2_old = t2_new
		
	def solveWithIntermediates(self):

		#Compute <phi0| H | phi0>
		Eref = self.computeEref()
		print "<phi0| H | phi0> = %.12f" % Eref
		#Initialize amplitudes
		self.initialize()
		
		self.updateIntermediates()
 
		Eold = self.ECCSD(self.t1_old,self.t2_old)
		print Eref+Eold
		
		#Compute 1 iteration
		t1_new = self.computeT1AmplitudesWithIntermediates(self.t1_old,self.t2_old)
		t2_new = self.computeT2AmplitudesWithIntermediates(self.t1_old,self.t2_old)
	
		Enew   = self.ECCSD(t1_new,t2_new)
		
		iters = 1
		print Eref+Enew, abs(Enew-Eold),iters
		self.t1_old = t1_new
		self.t2_old = t2_new

		
		while(abs(Enew-Eold) > self.precision):
			self.updateIntermediates()
			Eold   = Enew
			t1_new = self.computeT1AmplitudesWithIntermediates(self.t1_old,self.t2_old)
			t2_new = self.computeT2AmplitudesWithIntermediates(self.t1_old,self.t2_old)
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

	def computeT1AmplitudesWithIntermediates(self,t1,t2):

		N = self.holeStates
		L = self.BasisFunctions

		tsingles = np.zeros((L-N,N))

		for i in range(0,N):
			for a in range(N,L):
				
				tsingles[a-N,i] = self.F[a,i]

				tmp1 = 0
				for k in range(0,N):
					 tmp1 += t1[a-N,k]*self.F2[k,i]
				
				tsingles[a-N,i] -= tmp1	 
				

				for c in range(N,L):
					if(c != a):
						tsingles[a-N,i] += self.F[a,c]*t1[c-N,i]

				for k in range(0,N):
					for c in range(N,L):
						tsingles[a-N,i] += self.QRPS2(c,i,k,a)*t1[c-N,k]

				tmp2 = 0		
				for k in range(0,N):
					for l in range(0,N):
						for c in range(N,L):
							 tmp2 += t2[c-N,a-N,k,l]*self.W3[k,l,c-N,i] #POSSIBLE TYPO!!!!!!!!!!!!!!!!!!!!
				tsingles[a-N,i] -= 0.5*tmp2			 
				
				for k in range(0,N):
					for c in range(N,L):
						 tsingles[a-N,i] += t2[a-N,c-N,i,k]*self.F1[c-N,k]

				for k in range(0,N):
					for c in range(N,L):
						for d in range(N,L):
							tsingles[a-N,i] += 0.5*self.QRPS2(c,d,k,a)*t2[c-N,d-N,k,i]

				for k in range(0,N):
					for c in range(N,L):
						for d in range(N,L):
							tsingles[a-N,i] += self.QRPS2(c,d,k,a)*t1[c-N,k]*t1[d-N,i]

				Dia  = self.F[i,i] - self.F[a,a]
				tsingles[a-N,i] = tsingles[a-N,i]/Dia
		
		return tsingles

	def computeT2AmplitudesWithIntermediates(self,t1,t2):

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
								tdoubles[a-N,b-N,i,j] += 0.5*(t1[a-N,k]*t1[b-N,l] - t1[b-N,k]*t1[a-N,l] + t2[a-N,b-N,k,l])*self.W1[k,l,i,j]

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] -= t1[b-N,k]*self.W2[a-N,k,i,j]
							tdoubles[a-N,b-N,i,j] += t1[a-N,k]*self.W2[b-N,k,i,j]

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] += t2[a-N,b-N,j,k]*self.F2[k,i]
							tdoubles[a-N,b-N,i,j] -= t2[a-N,b-N,i,k]*self.F2[k,j]

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,d,a,b)*t2[c-N,d-N,i,j]

						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += t2[b-N,c-N,i,j]*self.F3[a-N,c-N]
							tdoubles[a-N,b-N,i,j] -= t2[a-N,c-N,i,j]*self.F3[b-N,c-N]

						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += self.QRPS2(c,j,a,b)*t1[c-N,i]
							tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,i,a,b)*t1[c-N,j]

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,j,a,b)*t1[c-N,i]
								tdoubles[a-N,b-N,i,j] -= 0.5*self.QRPS2(c,j,a,b)*t1[c-N,i]

						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += t2[b-N,c-N,j,k]*self.W4[a-N,k,i,c-N]
								tdoubles[a-N,b-N,i,j] -= t2[b-N,c-N,i,k]*self.W4[a-N,k,j,c-N]
								tdoubles[a-N,b-N,i,j] -= t2[a-N,c-N,j,k]*self.W4[b-N,k,i,c-N]
								tdoubles[a-N,b-N,i,j] += t2[a-N,c-N,i,k]*self.W4[b-N,k,j,c-N]

						Dabij = self.F[i,i] + self.F[j,j] - self.F[a,a] - self.F[b,b]
						tdoubles[a-N,b-N,i,j] /= Dabij
		return tdoubles

	def computeT2AmplitudesCCD(self,t2):

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
								tdoubles[a-N,b-N,i,j] += 0.5*t2[a-N,b-N,k,l]*self.W1[k,l,i,j]

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] += t2[a-N,b-N,j,k]*self.F2[k,i]
							tdoubles[a-N,b-N,i,j] -= t2[a-N,b-N,i,k]*self.F2[k,j]

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(c,d,a,b)*t2[c-N,d-N,i,j]

						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += t2[b-N,c-N,i,j]*self.F3[a-N,c-N]
							tdoubles[a-N,b-N,i,j] -= t2[a-N,c-N,i,j]*self.F3[b-N,c-N]

						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += t2[b-N,c-N,j,k]*self.W4[a-N,k,i,c-N]
								tdoubles[a-N,b-N,i,j] -= t2[b-N,c-N,i,k]*self.W4[a-N,k,j,c-N]
								tdoubles[a-N,b-N,i,j] -= t2[a-N,c-N,j,k]*self.W4[b-N,k,i,c-N]
								tdoubles[a-N,b-N,i,j] += t2[a-N,c-N,i,k]*self.W4[b-N,k,j,c-N]

						Dabij = self.F[i,i] + self.F[j,j] - self.F[a,a] - self.F[b,b]
						tdoubles[a-N,b-N,i,j] /= Dabij
		return tdoubles

	def computeW1(self,k,l,i,j):

		#Doublechecked: I

		N = self.holeStates
		L = self.BasisFunctions
		
		val = self.QRPS2(i,j,k,l)
		
		tmp = 0
		for c in range(N,L):
			val += self.QRPS2(c,j,k,l)*self.t1_old[c-N,i] - self.QRPS2(c,i,k,l)*self.t1_old[c-N,j]
			for d in range(N,L):
				tmp += self.QRPS2(c,d,k,l)*(self.t1_old[c-N,i]*self.t1_old[d-N,j] - self.t1_old[c-N,j]*self.t1_old[d-N,i] + self.t2_old[c-N,d-N,i,j])
		
		tmp *= 0.5
		val += tmp
		
		return val

	def computeW1CCD(self,k,l,i,j):

		N = self.holeStates
		L = self.BasisFunctions
		
		val = self.QRPS2(i,j,k,l)
		
		for c in range(N,L):
			for d in range(N,L):
				val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[c-N,d-N,i,j]
		return val

	def computeW2(self,a,k,i,j):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(i,j,a,k)
		tmp = 0
		for c in range(N,L):
			val += self.QRPS2(i,c,a,k)*self.t1_old[c-N,j] - self.QRPS2(j,c,a,k)*self.t1_old[c-N,i]
			for d in range(N,L):
				tmp += self.QRPS2(c,d,a,k)*(self.t2_old[c-N,d-N,i,j]+self.t1_old[c-N,i]*self.t1_old[d-N,j]-self.t1_old[c-N,j]*self.t1_old[d-N,i])
		tmp *= 0.5
		val += tmp
		return val

	def computeW2CCD(self,a,k,i,j):
		#Doublechecked: I 
		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(i,j,a,k)

		for c in range(N,L):
			for d in range(N,L):
				val += 0.5*self.QRPS2(c,d,a,k)*self.t2_old[c-N,d-N,i,j]
		return val
	
	def computeW3(self,k,l,c,i):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = 0
		val = self.QRPS2(c,i,k,l)

		for d in range(N,L):
			val += self.QRPS2(c,d,k,l)*self.t1_old[c-N,i]
		return val

	def computeW3CCD(self,k,l,c,i):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(c,i,k,l)
		
		return val

	def computeW4(self,a,k,i,c):

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(i,c,a,k)

		for d in range(N,L):
			val += self.QRPS2(d,c,a,k)*self.t1_old[d-N,i]

		tmp1 = 0	
		for l in range(0,N):
			tmp1 += self.t1_old[a-N,l]*self.W3[k,l,c-N,i]

		val -= tmp1	
		
		tmp2 = 0
		
		for l in range(0,N):
			for d in range(N,L):
				tmp2 += self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,i,l]
		
		tmp2 *= 0.5
		val  += tmp2

		return val

	def computeW4CCD(self,a,k,i,c):

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(i,c,a,k)
		for l in range(0,N):
			for d in range(N,L):
				val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,i,l]
		return val

	def computeW4manual(self,a,k,i,c):

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(i,c,a,k)

		for d in range(N,L):
			val += self.QRPS2(d,c,a,k)*self.t1_old[d-N,i]

		for l in range(0,N):
			val -= self.t1_old[a-N,l]*self.computeW3(k,l,c,i)

		for l in range(0,N):
			for d in range(N,L):
				val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,i,l]

		return val

	def computeF1(self,c,k):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.F[k,c]

		for l in range(0,N):
			for d in range(N,L):
				val += self.QRPS2(c,d,k,l)*self.t1_old[d-N,l]

		return val

	def computeF1CCD(self,c,k):

		#Doublechecked: I 
		N = self.holeStates
		L = self.BasisFunctions
		val = self.F[k,c]
		return val

	def computeF2(self,k,i):

		N = self.holeStates
		L = self.BasisFunctions

		val = (1.0-self.delta(k,i))*self.F[k,i]

		for c in range(N,L):
			val += self.t1_old[c-N,i]*self.computeF1(c,k)
			for l in range(0,N):
				val += self.QRPS2(i,c,k,l)*self.t1_old[c-N,l]
				for d in range(N,L):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[c-N,d-N,i,l]
		return val

	def computeF2CCD(self,k,i):

		N = self.holeStates
		L = self.BasisFunctions

		val = 0

		if(k != i):
			val = self.F[k,i]

		for c in range(N,L):
			for l in range(0,N):
				for d in range(N,L):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[c-N,d-N,i,l]
		return val

	def computeF3(self,a,c):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = 0
		if(c != a):
			val = self.F[a,c]

		for k in range(0,N):
			val -= self.t1_old[a-N,k]*self.F1[c-N,k]
			for d in range(N,L):
				val -= self.QRPS2(c,d,k,a)*self.t1_old[d-N,k]
				for l in range(0,N):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,k,l]

		return val

	def computeF3CCD(self,a,c):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = 0
		if(c != a):
			val = self.F[a,c]

		for k in range(0,N):
			for d in range(N,L):
				for l in range(0,N):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,k,l]

		return val

	def computeF3manual(self,a,c):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = 0
		if(c != a):
			val = self.F[a,c]

		for k in range(0,N):
			val -= self.t1_old[a-N,k]*self.computeF1(c,k)
			for d in range(N,L):
				val -= self.QRPS2(c,d,k,a)*self.t1_old[d-N,k]
				for l in range(0,N):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,k,l]
		return val


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
#inFile = open('HO_2d_10_nonzero.dat','r')
inFile = open('coulomb2.dat')
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
#ccsd_test.solve()
#ccsd_test.solveWithCCD()
ccsd_test.solveWithIntermediates()