import numpy as np

class CCSD:
	def __init__(self, N,L,w,oneBodyElements, precision=1e-8):
		
		self.holeStates   		= N
		self.BasisFunctions     = L
		self.particleStates     = L-N
		self.twoBodyIntegrals   = w
		self.oneBodyIntegrals 	= oneBodyElements
		self.F    				= np.zeros((L,L))
		self.t1_old 			= np.zeros((L-N,N))
		self.t2_old 			= np.zeros( (L-N,L-N,N,N) )
		self.precision   		= precision

		self.tau 				= np.zeros((L-N,L-N,N,N))
		self.tau2 				= np.zeros((L-N,L-N,N,N))
		self.W1 				= np.zeros((N,N,N,N))
		self.W2 				= np.zeros((L-N,L-N,L-N,L-N))
		self.W3 				= np.zeros((N,L-N,L-N,N))
		self.F1 				= np.zeros((N,L-N))
		self.F2 				= np.zeros((N,N))
		self.F3 				= np.zeros((L-N,L-N))

	def updateIntermediates(self):

		for a in range(N,L):
			for b in range(N,L):
				for i in range(0,N):
					for j in range(0,N):
						self.tau[a-N,b-N,i,j]  = self.computeTau(a,b,i,j)
						self.tau2[a-N,b-N,i,j] = self.computeTau2(a,b,i,j)

		for k in range(0,N):
			for l in range(0,N):
				for i in range(0,N):
					for j in range(0,N):
						self.W1[k,l,i,j] = self.computeW1(k,l,i,j)

		for a in range(N,L):
			for b in range(N,L):
				for e in range(N,L):
					for f in range(N,L):
						self.W2[a-N,b-N,e-N,f-N] = self.computeW2(a,b,e,f)

		for m in range(0,N):
			for b in range(N,L):
				for e in range(N,L):
					for j in range(0,N):
						self.W3[m,b-N,e-N,j] = self.computeW3(m,b,e,j)
		
		for c in range(N,L):
			for k in range(0,N):
				self.F1[k,c-N] = self.computeF1(k,c)

		for k in range(0,N):
			for i in range(0,N):
				self.F2[k,i] = self.computeF2(k,i)

		for a in range(N,L):
			for c in range(N,L):
				self.F3[a-N,c-N] = self.computeF3(a,c)
		
	def delta(self,i,j):
		if(i == j):
			return 1.0
		else:
			return 0

	def QRPS2(self,q,r,p,s):
		return self.twoBodyIntegrals[q,r,p,s]
	
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
		
		self.F = self.computeFockMatrix()

		N = self.holeStates
		L = self.BasisFunctions
		
		
		for i in range(0,N):
			for a in range(N,L):
				Dia  = self.F[i,i] - self.F[a,a]
				tia  = self.F[i,a]/Dia
				self.t1_old[a-N,i] = tia
		
		for i in range(0,N):
			for j in range(0,N):
				for a in range(N,L):
					for b in range(N,L):
						Dabij = self.F[i,i]+self.F[j,j] - self.F[a,a] - self.F[b,b]
						self.t2_old[a-N,b-N,i,j] = self.QRPS2(a,b,i,j)/Dabij

	def ECCSD(self,t1,t2):
		
		N = self.holeStates
		L = self.BasisFunctions

		ECCSD = 0.0
		
		for i in range(0,N):
			for a in range(N,L):
				ECCSD += self.F[i,a]*t1[a-N,i]


		tmp = 0		
		for i in range(0,N):
			for j in range(0,N):
				for a in range(N,L):
					for b in range(N,L):
						tmp += self.QRPS2(i,j,a,b)*(t2[a-N,b-N,i,j] + t1[a-N,i]*t1[b-N,j]-t1[b-N,i]*t1[a-N,j])

		return ECCSD+0.25*tmp

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
					
	def solveWithIntermediates(self):

		#Compute <phi0| H | phi0>
		Eref = self.computeEref()
		print "<phi0| H | phi0> = %g" % Eref
		#Initialize amplitudes
		self.initialize()

		self.updateIntermediates()
 		
		Eold = self.ECCSD(self.t1_old,self.t2_old)
		print Eref+Eold, Eold
		
		#Compute 1 iteration
		t1_new = self.computeT1SG(self.t1_old,self.t2_old)
		t2_new = self.computeT2SG(self.t1_old,self.t2_old)
		
		Enew   = self.ECCSD(t1_new,t2_new)
		
		iters = 1
		print Eref+Enew, abs(Enew-Eold),iters
		self.t1_old = t1_new
		self.t2_old = t2_new

		while(abs(Enew-Eold) > self.precision):
			self.updateIntermediates()
			Eold   = Enew
			t1_new = self.computeT1SG(self.t1_old,self.t2_old)
			t2_new = self.computeT2SG(self.t1_old,self.t2_old)
			Enew   = self.ECCSD(t1_new,t2_new)
			
			iters += 1
			print Eref+Enew, abs(Enew-Eold), iters
			self.t1_old = t1_new
			self.t2_old = t2_new

	def computeT1SG(self,t1,t2):

		N = self.holeStates
		L = self.BasisFunctions

		tsingles = np.zeros((L-N,N))

		for i in range(0,N):
			for a in range(N,L):
				tsingles[a-N,i] = self.computeZ1(a,i)/(self.F[i,i]-self.F[a,a])

		return tsingles

	def computeT2SG(self,t1,t2):

		N = self.holeStates
		L = self.BasisFunctions

		tdoubles = np.zeros( (L-N,L-N,N,N) )

		for a in range(N,L):
			for b in range(N,L):
				for i in range(0,N):
					for j in range(0,N):
						Dabij = self.F[i,i] + self.F[j,j] - self.F[a,a] - self.F[b,b]
						tdoubles[a-N,b-N,i,j] = self.computeZ2(a,b,i,j)/Dabij

		return tdoubles

	def computeZ1(self,a,i):
		
		val = self.F[i,a]

		for e in range(N,L):
			val += self.t1_old[e-N,i]*self.F3[a-N,e-N]

		for m in range(0,N):
			val -= self.t1_old[a-N,m]*self.F2[m,i]
			for e in range(N,L):
				val += self.t2_old[a-N,e-N,i,m]*self.F1[m,e-N]

		for n in range(0,N):
			for f in range(N,L):
				val -= self.t1_old[f-N,n]*self.QRPS2(n,a,i,f)

		for m in range(0,N):
			for e in range(N,L):
				for f in range(N,L):
					val -= 0.5*self.t2_old[e-N,f-N,i,m]*self.QRPS2(m,a,e,f)

		for m in range(0,N):
			for e in range(N,L):
				for n in range(0,N):
					val -= 0.5*self.t2_old[a-N,e-N,m,n]*self.QRPS2(n,m,e,i)
		
		return val

	def computeZ2(self,a,b,i,j):

		val = self.QRPS2(i,j,a,b)

		for e in range(N,L):
			val += self.t2_old[a-N,e-N,i,j]*self.F3[b-N,e-N] - self.t2_old[b-N,e-N,i,j]*self.F3[a-N,e-N]
			for m in range(0,N):
				val -= 0.5*self.t2_old[a-N,e-N,i,j]*self.t1_old[b-N,m]*self.F1[m,e-N]
				val += 0.5*self.t2_old[b-N,e-N,i,j]*self.t1_old[a-N,m]*self.F1[m,e-N]


		sum2 = 0
		for m in range(0,N):
			sum2 += self.t2_old[a-N,b-N,i,m]*self.F2[m,j] - self.t2_old[a-N,b-N,j,m]*self.F2[m,i]
			for e in range(N,L):
				sum2 += 0.5*self.t1_old[e-N,j]*self.t2_old[a-N,b-N,i,m]*self.F1[m,e-N]
				sum2 -= 0.5*self.t1_old[e-N,i]*self.t2_old[a-N,b-N,j,m]*self.F1[m,e-N] 
		val -= sum2

		for m in range(0,N):
			for n in range(0,N):
				val += 0.5*self.tau[a-N,b-N,m,n]*self.W1[m,n,i,j]

		for e in range(N,L):
			for f in range(N,L):
				val += 0.5*self.tau[e-N,f-N,i,j]*self.W2[a-N,b-N,e-N,f-N]

		for m in range(0,N):
			for e in range(N,L):
				val += (self.t2_old[a-N,e-N,i,m]*self.W3[m,b-N,e-N,j] - self.t1_old[e-N,i]*self.t1_old[a-N,m]*self.QRPS2(m,b,e,j))
				val -= (self.t2_old[b-N,e-N,i,m]*self.W3[m,a-N,e-N,j] - self.t1_old[e-N,i]*self.t1_old[b-N,m]*self.QRPS2(m,a,e,j))
				val -= (self.t2_old[a-N,e-N,j,m]*self.W3[m,b-N,e-N,i] - self.t1_old[e-N,j]*self.t1_old[a-N,m]*self.QRPS2(m,b,e,i))
				val += (self.t2_old[b-N,e-N,j,m]*self.W3[m,a-N,e-N,i] - self.t1_old[e-N,j]*self.t1_old[b-N,m]*self.QRPS2(m,a,e,i))

		for e in range(N,L):
			val += self.t1_old[e-N,i]*self.QRPS2(a,b,e,j) - self.t1_old[e-N,j]*self.QRPS2(a,b,e,i)

		for m in range(0,N):
			val -= self.t1_old[a-N,m]*self.QRPS2(m,b,i,j) - self.t1_old[b-N,m]*self.QRPS2(m,a,i,j)
		return val
			
	def computeT1Amplitudes(self,t1,t2):
		
		N = self.holeStates
		L = self.BasisFunctions

		tsingles = np.zeros((L-N,N))

		for i in range(0,N):
			for a in range(N,L):
				
				tsingles[a-N,i] = self.F[a,i]

				#Crawford and Schaefer

				for b in range(N,L):
					tsingles[a-N,i] += (1.0-self.delta(b,a))*self.F[a,b]*t1[b-N,i]
				for j in range(0,N):
					tsingles[a-N,i] -= (1.0-self.delta(i,j))*self.F[j,i]*t1[a-N,j]

				"""
				for j in range(0,N):
					for b in range(N,L):
						tsingles[a-N,i] += self.QRPS2(i,a,b,j)*t1[b-N,i]

						tsingles[a-N,i] += self.F[i,a]*t2[b-N,a-N,j,i]

						tsingles[a-N,i] -= self.F[i,a]*t1[a-N,j]*t1[b-N,i]

				for j in range(0,N):
					for k in range(0,N):
						for b in range(N,L):
							tsingles[a-N,i] -= 0.5*self.QRPS2(i,j,k,a)*t2[b-N,a-N,i,j]

							tsingles[a-N,i] += self.QRPS2(i,j,a,k)*t1[b-N,i]*t1[a-N,j]

				for j in range(0,N):
					for b in range(N,L):
						for c in range(N,L):
							tsingles[a-N,i] -= self.QRPS2(i,a,b,c)*t1[b-N,j]*t1[c-N,i]

							tsingles[a-N,i] += 0.5*self.QRPS2(a,i,b,c)*t2[b-N,c-N,j,i]

				for j in range(0,N):
					for k in range(0,N):
						for b in range(N,L):
							for c in range(N,L):
								tsingles[a-N,i] += self.QRPS2(i,j,a,b)*t1[a-N,i]*t2[b-N,c-N,j,k]

								tsingles[a-N,i] += 0.5*self.QRPS2(i,j,a,b)*t1[a-N,k]*t2[b-N,c-N,i,j]

								tsingles[a-N,i] += 0.5*self.QRPS2(i,j,a,b)*t1[c-N,i]*t2[a-N,b-N,j,k]

								tsingles[a-N,i] += self.QRPS2(i,j,a,b)*t1[a-N,k]*t1[b-N,i]*t1[c-N,j]
				"""				
				
				for k in range(0,N):
					for c in range(N,L):
						tsingles[a-N,i] += self.QRPS2(k,a,c,i)*t1[c-N,k]

						tsingles[a-N,i] += self.F[k,c]*t2[a-N,c-N,i,k]

						tsingles[a-N,i] -= self.F[k,c]*t1[c-N,i]*t1[a-N,k]

				for k in range(0,N):
					for c in range(N,L):
						for d in range(N,L):
							kacd = self.QRPS2(k,a,c,d)
							tsingles[a-N,i] += 0.5*kacd*t2[c-N,d-N,k,i]

							tsingles[a-N,i] -= kacd*t1[c-N,k]*t1[d-N,i]

				for k in range(0,N):
					for l in range(0,N):
						for c in range(N,L):

							klci = self.QRPS2(k,l,c,i)
							
							tsingles[a-N,i] -= 0.5*klci*t2[c-N,a-N,k,l]

							tsingles[a-N,i] -= klci*t1[c-N,k]*t1[a-N,l]

				for k in range(0,N):
					for l in range(0,N):
						for c in range(N,L):
							for d in range(N,L):

								klcd = self.QRPS2(k,l,c,d)

								tsingles[a-N,i] -= klcd*t1[c-N,k]*t1[d-N,i]*t1[a-N,l]

								tsingles[a-N,i] += klcd*t1[c-N,k]*t2[d-N,a-N,l,i]

								tsingles[a-N,i] -= 0.5*klcd*t2[c-N,d-N,k,i]*t1[a-N,l]

								tsingles[a-N,i] -= 0.5*klcd*t2[c-N,a-N,k,l]*t1[d-N,i]

				
				"""
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
				"""				
				
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
						
						tdoubles[a-N,b-N,i,j] = self.QRPS2(a,b,i,j)
						###########################################Crawford and Schaefer
						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] -= (1.0-self.delta(k,j))*self.F[k,j]*t2[a-N,b-N,i,k] - (1-self.delta(k,i))*self.F[k,i]*t2[a-N,b-N,j,k]

						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += (1.0-self.delta(b,c))*self.F[b,c]*t2[a-N,c-N,i,j] - (1-self.delta(a,c))*self.F[a,c]*t2[b-N,c-N,i,j]
						


						for k in range(0,N):
							for l in range(0,N):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(k,l,i,j)*t2[a-N,b-N,k,l]

								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(k,l,i,j)*(t1[a-N,k]*t1[b-N,l] - t1[b-N,k]*t1[a-N,l])

						for c in range(N,L):
							for d in range(N,L):
								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(a,b,c,d)*t2[c-N,d-N,i,j]

								tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(a,b,c,d)*(t1[c-N,i]*t1[d-N,j] - t1[c-N,j]*t1[d-N,i])

						################## New stuff		

						for k in range(0,N):
							for c in range(N,L):
								
								tdoubles[a-N,b-N,i,j] += self.QRPS2(k,b,c,j)*t2[a-N,c-N,i,k] - self.QRPS2(k,a,c,j)*t2[b-N,c-N,i,k]
								tdoubles[a-N,b-N,i,j] -= self.QRPS2(k,b,c,i)*t2[a-N,c-N,j,k] - self.QRPS2(k,a,c,i)*t2[b-N,c-N,j,k]

								tdoubles[a-N,b-N,i,j] -= self.QRPS2(k,b,i,c)*t1[a-N,k]*t1[c-N,j] - self.QRPS2(k,a,i,c)*t1[b-N,k]*t1[c-N,j]
								tdoubles[a-N,b-N,i,j] += self.QRPS2(k,b,j,c)*t1[a-N,k]*t1[c-N,i] - self.QRPS2(k,a,j,c)*t1[b-N,k]*t1[c-N,i]

								tdoubles[a-N,b-N,i,j] += self.F[k,c]*(t1[a-N,k]*t2[b-N,c-N,i,j] - t1[b-N,k]*t2[a-N,c-N,i,j])

								tdoubles[a-N,b-N,i,j] += self.F[k,c]*(t1[c-N,i]*t2[a-N,b-N,j,k] - t1[c-N,j]*t2[a-N,b-N,i,k])
	

						for c in range(N,L):
							tdoubles[a-N,b-N,i,j] += self.QRPS2(a,b,c,j)*t1[c-N,i] - self.QRPS2(a,b,c,i)*t1[c-N,j]

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] -= self.QRPS2(k,b,i,j)*t1[a-N,k] - self.QRPS2(k,a,i,j)*t1[b-N,k]

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										
										klcd   = self.QRPS2(k,l,c,d)
										klcd_2 = (1.0/2.0)*klcd
										klcd_4 = (1.0/4.0)*klcd

										tdoubles[a-N,b-N,i,j] += klcd_2*(t2[a-N,c-N,i,k]*t2[d-N,b-N,l,j] - t2[b-N,c-N,i,k]*t2[d-N,a-N,l,j])
										tdoubles[a-N,b-N,i,j] -= klcd_2*(t2[a-N,c-N,j,k]*t2[d-N,b-N,l,i] - t2[b-N,c-N,j,k]*t2[d-N,a-N,l,i])

										tdoubles[a-N,b-N,i,j] += klcd_4*t2[c-N,d-N,i,j]*t2[a-N,b-N,k,l]

										tdoubles[a-N,b-N,i,j] -= klcd_2*(t2[a-N,c-N,i,j]*t2[b-N,d-N,k,l] - t2[b-N,c-N,i,j]*t2[a-N,d-N,k,l])

										tdoubles[a-N,b-N,i,j] -= klcd_2*(t2[a-N,b-N,i,k]*t2[c-N,d-N,j,l] - t2[a-N,b-N,j,k]*t2[c-N,d-N,i,l])							

										tdoubles[a-N,b-N,i,j] -= klcd*(t1[c-N,k]*t1[d-N,i]*t2[a-N,b-N,l,j] - t1[c-N,k]*t1[d-N,j]*t2[a-N,b-N,l,i])

										tdoubles[a-N,b-N,i,j] -= klcd*(t1[c-N,k]*t1[a-N,l]*t2[d-N,b-N,i,j] - t1[c-N,k]*t1[b-N,l]*t2[d-N,a-N,i,j])

										tdoubles[a-N,b-N,i,j] += klcd_4*(t1[c-N,i]*t1[d-N,j]*t2[a-N,b-N,k,l] - t1[c-N,j]*t1[d-N,i]*t2[a-N,b-N,k,l])

										tdoubles[a-N,b-N,i,j] += klcd_4*(t1[a-N,k]*t1[b-N,l]*t2[c-N,d-N,i,j] - t1[b-N,k]*t1[a-N,l]*t2[c-N,d-N,i,j])

										tdoubles[a-N,b-N,i,j] += klcd*(t1[c-N,i]*t1[b-N,l]*t2[a-N,d-N,k,j] - t1[c-N,i]*t1[a-N,l]*t2[b-N,d-N,k,j])
										tdoubles[a-N,b-N,i,j] -= klcd*(t1[c-N,j]*t1[b-N,l]*t2[a-N,d-N,k,i] - t1[c-N,j]*t1[a-N,l]*t2[b-N,d-N,k,i])

										tdoubles[a-N,b-N,i,j] += klcd_4*(t1[c-N,i]*t1[a-N,k]*t1[d-N,j]*t1[b-N,l] - t1[c-N,i]*t1[b-N,k]*t1[d-N,j]*t1[a-N,l])
										tdoubles[a-N,b-N,i,j] -= klcd_4*(t1[c-N,j]*t1[a-N,k]*t1[d-N,i]*t1[b-N,l] - t1[c-N,j]*t1[b-N,k]*t1[d-N,i]*t1[a-N,l])
						
								
						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									tdoubles[a-N,b-N,i,j] -= self.QRPS2(k,l,c,i)*t1[c-N,k]*t2[a-N,b-N,l,j] - self.QRPS2(k,l,c,j)*t1[c-N,k]*t2[a-N,b-N,l,i]

									tdoubles[a-N,b-N,i,j] += self.QRPS2(k,l,i,c)*t1[a-N,l]*t2[b-N,c-N,j,k] - self.QRPS2(k,l,i,c)*t1[b-N,l]*t2[a-N,c-N,j,k]
									tdoubles[a-N,b-N,i,j] -= self.QRPS2(k,l,j,c)*t1[a-N,l]*t2[b-N,c-N,i,k] - self.QRPS2(k,l,j,c)*t1[b-N,l]*t2[a-N,c-N,i,k]

									tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(k,l,c,j)*t1[c-N,i]*t2[a-N,b-N,k,l] - 0.5*self.QRPS2(k,l,c,i)*t1[c-N,j]*t2[a-N,b-N,k,l]

									tdoubles[a-N,b-N,i,j] += 0.5*self.QRPS2(k,l,c,j)*t1[c-N,i]*t1[a-N,k]*t1[b-N,l] - 0.5*self.QRPS2(k,l,c,j)*t1[c-N,i]*t1[b-N,k]*t1[a-N,l]
									tdoubles[a-N,b-N,i,j] -= 0.5*self.QRPS2(k,l,c,i)*t1[c-N,j]*t1[a-N,k]*t1[b-N,l] - 0.5*self.QRPS2(k,l,c,i)*t1[c-N,j]*t1[b-N,k]*t1[a-N,l]

						for k in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									tdoubles[a-N,b-N,i,j] += self.QRPS2(k,a,c,d)*t1[c-N,k]*t2[d-N,b-N,i,j] - self.QRPS2(k,b,c,d)*t1[c-N,k]*t2[d-N,a-N,i,j]

									tdoubles[a-N,b-N,i,j] += self.QRPS2(a,k,d,c)*t1[d-N,i]*t2[b-N,c-N,j,k] - self.QRPS2(b,k,d,c)*t1[d-N,i]*t2[a-N,c-N,j,k]
									tdoubles[a-N,b-N,i,j] -= self.QRPS2(a,k,d,c)*t1[d-N,j]*t2[b-N,c-N,i,k] - self.QRPS2(b,k,d,c)*t1[d-N,j]*t2[a-N,c-N,i,k]

									kbcd = self.QRPS2(k,b,c,d); kacd = self.QRPS2(k,a,c,d)

									tdoubles[a-N,b-N,i,j] -= 0.5*(kbcd*t1[a-N,k]*t2[c-N,d-N,i,j] - kacd*t1[b-N,k]*t2[c-N,d-N,i,j])

									tdoubles[a-N,b-N,i,j] -= 0.5*(kbcd*t1[c-N,i]*t1[a-N,k]*t1[d-N,j] - kacd*t1[c-N,i]*t1[b-N,k]*t1[d-N,j])
									tdoubles[a-N,b-N,i,j] += 0.5*(kbcd*t1[c-N,j]*t1[a-N,k]*t1[d-N,i] - kacd*t1[c-N,j]*t1[b-N,k]*t1[d-N,i])

						################################ Old stuff
						"""
						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										tdoubles[a-N,b-N,i,j] += 0.25*self.QRPS2(c,d,k,l)*t2[a-N,b-N,k,l]

						

						#P(ab) -terms
						for k in range(0,N):
							for c in range(N,L):
								tdoubles[a-N,b-N,i,j] += self.F[k,c]*(t1[a-N,k]*t2[b-N,c-N,i,j] - t1[b-N,k]*t2[a-N,c-N,i,j])

						for k in range(0,N):
							tdoubles[a-N,b-N,i,j] -= self.QRPS2(i,j,k,b)*t1[a-N,k] - self.QRPS2(i,j,k,a)*t1[b-N,k]

						for k in range(0,N):
							for l in range(0,N):
								tdoubles[a-N,b-N,i,j] += 0.5*(self.QRPS2(i,j,k,l)*(t1[a-N,k]*t1[b-N,l] - t1[b-N,k]*t1[a-N,l]))

						for k in range(0,N):
							for c in range(N,L):
								for d in range(N,L):
									  tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,a)*(t1[c-N,k]*t2[d-N,b-N,i,j]  + 0.5*t1[b-N,k]*t2[c-N,d-N,i,j])
									  tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,b)*(t1[c-N,k]*t2[d-N,a-N,i,j]  + 0.5*t1[a-N,k]*t2[c-N,d-N,i,j])

						for k in range(0,N):
							for l in range(0,N):
								for c in range(N,L):
									for d in range(N,L):
										tdoubles[a-N,b-N,i,j] += self.QRPS2(c,d,k,l)*(0.25*t1[a-N,k]*t1[b-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[a-N,l]*t2[d-N,b-N,i,j] - 0.5*t2[a-N,c-N,i,j]*t2[b-N,d-N,k,l])
										tdoubles[a-N,b-N,i,j] -= self.QRPS2(c,d,k,l)*(0.25*t1[b-N,k]*t1[a-N,l]*t2[c-N,d-N,i,j] - t1[c-N,k]*t1[b-N,l]*t2[d-N,a-N,i,j] - 0.5*t2[b-N,c-N,i,j]*t2[a-N,d-N,k,l])

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
						"""				
						Dabij = self.F[i,i] + self.F[j,j] - self.F[a,a] - self.F[b,b]
						tdoubles[a-N,b-N,i,j] /= Dabij
		return tdoubles

	def computeTau(self,a,b,i,j):

		val = self.t2_old[a-N,b-N,i,j] + self.t1_old[a-N,i]*self.t1_old[b-N,j] - self.t1_old[b-N,i]*self.t1_old[a-N,j]
		return val

	def computeTau2(self,a,b,i,j):
		val = self.t2_old[a-N,b-N,i,j] + 0.5*(self.t1_old[a-N,i]*self.t1_old[b-N,j] - self.t1_old[b-N,i]*self.t1_old[a-N,j])
		return val

	def computeW1(self,m,n,i,j):

		#Doublechecked: I

		N = self.holeStates
		L = self.BasisFunctions
		
		val = self.QRPS2(m,n,i,j)

		for e in range(N,L):
			val += self.QRPS2(m,n,i,e)*self.t1_old[e-N,j] - self.QRPS2(m,n,j,e)*self.t1_old[e-N,i]
			for f in range(N,L):
				val += 0.25*self.QRPS2(m,n,e,f)*self.tau[e-N,f-N,i,j]

		return val

	def computeW2(self,a,b,e,f):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(a,b,e,f)

		for m in range(0,N):
			val -= self.t1_old[b-N,m]*self.QRPS2(a,m,e,f) - self.t1_old[a-N,m]*self.QRPS2(b,m,e,f)
			for n in range(0,N):
				val += 0.25*self.tau[a-N,b-N,m,n]*self.QRPS2(m,n,e,f)
		
		return val

	def computeW3(self,m,b,e,j):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.QRPS2(m,b,e,j)

		for f in range(N,L):
			val += self.QRPS2(m,b,e,f)*self.t1_old[f-N,j]
		
		for n in range(0,N):
			val -= self.t1_old[b-N,n]*self.QRPS2(m,n,e,j)

		for n in range(0,N):
			for f in range(N,L):
				val -= (0.5*self.t2_old[f-N,b-N,j,n]+self.t1_old[f-N,j]*self.t1_old[b-N,n])*self.QRPS2(m,n,e,f)

		return val

	def computeF1(self,m,e):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = self.F[m,e]

		for n in range(0,N):
			for f in range(N,L):
				val += self.QRPS2(m,n,e,f)*self.t1_old[f-N,n]

		return val

	def computeF2(self,m,i):

		N = self.holeStates
		L = self.BasisFunctions

		val = (1.0-self.delta(m,i))*self.F[m,i]
		
		for e in range(N,L):
			val += 0.5*self.F1[m,e-N]*self.t1_old[e-N,i]

		for e in range(N,L):
			for n in range(0,N):
				val += self.QRPS2(m,n,i,e)*self.t1_old[e-N,n]

		for n in range(0,N):
			for e in range(N,L):
				for f in range(N,L):
					val += 0.5*self.QRPS2(m,n,e,f)*self.tau2[e-N,f-N,i,n]

		return val

	def computeF3(self,a,e):

		#Doublechecked: I 

		N = self.holeStates
		L = self.BasisFunctions

		val = (1.0-self.delta(a,e))*self.F[a,e]

		for m in range(0,N):
			val -= 0.5*self.F[e-N,m]*self.t1_old[a-N,m]

		for m in range(0,N):
			for n in range(0,N):
				for f in range(N,L):
					val -= 0.5*self.QRPS2(m,n,e,f)*self.tau2[a-N,f-N,m,n]

		for m in range(0,N):
			for f in range(N,L):
				val += self.QRPS2(m,a,f,e)*self.t1_old[f-N,m]

		"""
		for k in range(0,N):
			val -= self.t1_old[a-N,k]*self.F1[c-N,k]
			for d in range(N,L):
				val -= self.QRPS2(c,d,k,a)*self.t1_old[d-N,k]
				for l in range(0,N):
					val += 0.5*self.QRPS2(c,d,k,l)*self.t2_old[a-N,d-N,k,l]
		"""
		return val

def delta(i,j):
	if(i == j):
		return 1
	else:
		return 0
def map_index_SPenergy(x):
	
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
def map_index2(x):
	
	if(x == 0 or x == 1):
		return 1

	elif(x == 2 or x == 3):
		return 2

	elif (x == 4 or x == 5):
		return 3

	elif (x == 6 or x == 7):
		return 4

	elif (x == 8 or x == 9):
		return 5

	elif (x == 10 or x == 11):
		return 6

	elif (x == 12 or x == 13):
		return 7

	elif (x == 14 or x == 15):
		return 8

	elif (x == 16 or x == 17):
		return 9

	elif (x == 18 or x == 19):
		return 10
def SingleParticleEnergy(x):
	if (x == 0):
		return 1.0

	elif (x == 1 or x == 2):
		return 2.0

	elif (x == 3 or x == 4 or x == 5):
		return 3.0

	elif (x == 6 or x == 7 or x == 8 or x == 9):
		return 4.0
def spin(x):
	if(   (x == 0) | (x == 2) | (x == 4) | (x == 6) | (x == 8) | (x == 10) | (x == 12) | (x == 14) | (x == 16) | (x == 18) ):
		return -1
	elif( (x == 1) | (x == 3) | (x == 5) | (x == 7) | (x == 9) | (x == 11) | (x == 13) | (x == 15) | (x == 17) | (x == 19) ):
		return 1

def QRPS(twoBodyIntegrals,q,r,p,s):
	
		key1 = str(map_index(q)) +str(map_index(r)) + str(map_index(p)) + str(map_index(s))
		key2 = str(map_index(q)) +str(map_index(r)) + str(map_index(s)) + str(map_index(p))

		d_qp = delta(spin(q),spin(p))
		d_rs = delta(spin(r),spin(s))
		d_qs = delta(spin(q),spin(s))
		d_rp = delta(spin(r),spin(p))

		if((twoBodyIntegrals.has_key(key1) == True) and (twoBodyIntegrals.has_key(key2) == True)):
			qrps = d_qp*d_rs*twoBodyIntegrals[key1]-d_qs*d_rp*twoBodyIntegrals[key2]
		elif ((twoBodyIntegrals.has_key(key1) == True) and (twoBodyIntegrals.has_key(key2) == False)):
			qrps = d_qp*d_rs*twoBodyIntegrals[key1]
		elif ((twoBodyIntegrals.has_key(key1) == False) and (twoBodyIntegrals.has_key(key2) == True)):
			qrps = -d_qs*d_rp*twoBodyIntegrals[key2]
		else:
			qrps = 0
			
		return qrps

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
w2 = np.zeros((L,L,L,L))

for i in range(0,L):
	for j in range(0,L):
		for k in range(0,L):
			for l in range(0,L):
				w2[i,j,k,l] = QRPS(w,i,j,k,l)

w = {}

prec = 1e-8

oneBodyElements = np.zeros(L)

for i in range(0,L):
	oneBodyElements[i] = SingleParticleEnergy(map_index_SPenergy(i))

ccsd_test = CCSD(N,L,w2,oneBodyElements)
#ccsd_test.solve()
#ccsd_test.solveWithCCD()
ccsd_test.solveWithIntermediates()

#####
#Brute force with t1 set to zero at all iterations:

#3.253314137
#3.15968739638
#3.14699374557 0.0126936508082 1
#3.14349585212 0.00349789345545 2
#3.14236761366 0.00112823845747 3
#3.14200235656 0.000365257095479 4
#3.14188369589 0.000118660672445 5
#3.14184504961 3.86462809903e-05 6
#3.14183244118 1.26084298357e-05 7
#3.14182832286 4.11831901667e-06 8
#3.14182697664 1.34621628466e-06 9
#3.14182653636 4.40282512532e-07 10
#3.14182639232 1.44043596231e-07 11
#3.14182634518 4.71359577064e-08 12
#3.14182632976 1.54267210117e-08 13
#3.14182632471 5.0493580428e-09 14

#real	3m59.862s
#user	3m58.624s
#sys	0m0.992s

#With Intermediates with t1 set to zero at all iterations:

#3.253314137
#3.15968739638
#3.14699374557 0.0126936508082 1
#3.14349585212 0.00349789345545 2
#3.14236761366 0.00112823845747 3
#3.14200235656 0.000365257095479 4
#3.14188369589 0.000118660672445 5
#3.14184504961 3.86462809904e-05 6
#3.14183244118 1.26084298356e-05 7
#3.14182832286 4.11831901662e-06 8
#3.14182697664 1.34621628452e-06 9
#3.14182653636 4.40282512643e-07 10
#3.14182639232 1.44043596287e-07 11
#3.14182634518 4.71359575954e-08 12
#3.14182632976 1.54267209701e-08 13
#3.14182632471 5.04935813994e-09 14