from numpy import *

def energy(state):
	if(1<=state<=2):
		return 1.0
	elif (3<=state<=6):
		return 2.0
	elif (7<=state<=12):
		return 3.0

def delta(i,j):
	if(i == j):
		return 1
	else:
		return 0
def mapping(state):
	if((state == 1) | (state==2)):
		return "1" #(n=0, m=0)
	elif((state == 3) | (state == 4)):
		return "2" #(n=0,m=-1)
	elif ((state == 5) | (state == 6)):
		return "3" #(n=0,m=1)
	elif ((state == 7) | (state == 8)):
		return "4" #(n=0,m=-2)
	elif ((state==9) | (state == 10)):
		return "5" #(n=1, m=0)
	elif ((state == 11) | (state == 12)):
		return "6"	#(n=0,m=2)

def spin(state):
	if((state%2 == 0)):
		return -2 #Spin down if the state is an even number
	else:
		return 2 #Spin up if the state is an odd number	
		
def expref(N,w):
	#<ref|H|ref>
	E0 = 0.0
	for i in range(1,N+1):
		E0+=energy(i)
		for j in range(1,N+1):
			key1 = mapping(i) + mapping(j) + mapping(i) + mapping(j)
			key2 = mapping(i) + mapping(j) + mapping(j) + mapping(i)
			E0+= 0.5*(delta(spin(i),spin(i))*delta(spin(j),spin(j))*w[key1] - delta(spin(i),spin(j))*delta(spin(j),spin(i))*w[key2])
	return E0
			
def exp_ref_refia(i,a,N,w):
	#<ref|H|ref_i^a>
	val = 0.0
	for j in range(1,N+1):
		key1 = mapping(i) + mapping(j) + mapping(a) + mapping(j)
		key2 = mapping(i) + mapping(j) + mapping(j) + mapping(a)
		val+= delta(spin(i),spin(a))*delta(spin(j),spin(j))*w[key1] - delta(spin(i),spin(j))*delta(spin(j),spin(a))*w[key2]
		
	return val

inFile = open('coulomb.dat','r')

w = {}

for line in inFile:
	tmp = line.split()
	key = tmp[0] + tmp[1] + tmp[2] + tmp[3]
	val = float(tmp[4])
	w[key] = val
	
	
inFile.close()

N = 2

exp_ref = expref(N,w)
exp_ref_19 = exp_ref_refia(1,9,N,w)
exp_ref_210 = exp_ref_refia(2,10,N,w)

print exp_ref




		
	


