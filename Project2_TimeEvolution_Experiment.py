from numpy import *
from scipy.linalg import expm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time

def V(x):
    return x**2
def Coulomb(w,p):
    return (w**2)*(p**2) #+1.0/p

pmax  = 20.0
pmin  = 0.0
Nstep = 200
h = (pmax-pmin)/float(Nstep)

A = zeros((Nstep,Nstep))
p = zeros(Nstep)
p[0] = pmin
p[Nstep-1] = pmax

for i in range(1,Nstep-1):
    p[i] = i*h
    A[i,i] = 2.0/(h**2) + Coulomb(1.0,p[i])
    
    if i < Nstep-2:
        A[i,i+1] = -1.0/(h**2)
        A[i+1,i] = -1.0/(h**2)
        
start_time = time.time()     
eig_vals, eig_vecs = linalg.eigh(A)
end_time = time.time()

print "Time: %.2f seconds" % (end_time-start_time)

eig_vals = sort(eig_vals)
eig_vecs = eig_vecs[:, eig_vals.argsort()]

#Time propagation
tfinal = 1.5
dt = 0.00001
Nt = int(tfinal/dt)
E0 = 1.0 #Strength of electric field
w = (eig_vals[3]-eig_vals[2])

psi_old = eig_vecs[:,2]
plt.figure(1)
plt.axis([0,pmax,0,0.15])
plt.plot(p,abs(psi_old**2))
plt.savefig("init.png")
norm    = sum(abs(psi_old)**2)
print "Enew=%g, norm=%g, t=%g" % (dot(conj(psi_old),dot(A,psi_old)), sum(abs(psi_old)**2),0)
n = 1

E_vec = [dot(conj(psi_old),dot(A,psi_old))]
t_vec = [0]

while ( (n <= Nt) and (norm < 1.0005) ):
    t  = n*dt
    t_vec.append(t)
    Ht = zeros((Nstep,Nstep))
    H1 = zeros((Nstep,Nstep))
    
    for i in range(1,Nstep-1):
        H1[i,i] = p[i]*E0*sin(w*t)
    
    Ht = A+H1
    
    Htilde  = eye(Nstep) - 1j*dt*Ht

    psi_new = dot(Htilde,psi_old) 
    E_vec.append(dot(conj(psi_new),dot(Ht,psi_new)).real)
    norm    = sum(abs(psi_new)**2)
    
    if(n % 10000 == 0):
        print t, norm
        plt.figure(n)
        plt.axis([0,pmax,0,0.15])
        plt.plot(p,abs(psi_new**2))
        plt.savefig("test%d.png" % n)
        #print "Enew=%g, norm=%g, t=%g" % (dot(conj(psi_new),dot(Ht,psi_new).real), sum(abs(psi_new)**2),t)
    
    psi_old = psi_new
    n += 1

t_vec = array(t_vec)
E_vec = array(E_vec)

plt.figure(2)
plt.plot(t_vec,E_vec)
plt.savefig("Energies.png")
#plt.show()
"""
energies_t = zeros(Nt+1)
energies_t[0] = eig_vals[2]

fig = plt.figure(1)
ax = fig.add_subplot(111) # Add axes
ax.grid()

line, = ax.plot(p, abs(psi_old)**2,'-r') # Fetch the line object
time_template = 'time = %.2fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
ax.set_xlabel('$p$') # X label of axes
ax.set_ylabel('$|\Psi(x, t)|^2$') # Y label of axes
ax.set_xlim(pmin, pmax)
ax.set_ylim(0,1)

H0 = zeros((Nstep,Nstep))


for i in range(1,Nstep-1):
    H0[i,i] = eig_vals[i+1]

print "Enew=%g" % (dot(conj(psi_old),dot(A,psi_old)))

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def animate(n):
    global energies_t
    global psi_old
    #global ck_old
    t  = n*dt
    Ht = zeros((Nstep,Nstep))
    H1 = zeros((Nstep,Nstep))
    
    for i in range(1,Nstep-1):
        H1[i,i] = E0*sin(w*t)
    
    Ht = A+H1
    
    Htilde  = eye(Nstep) - 1j*dt*Ht

    psi_new = dot(Htilde,psi_old) 

    print "Enew=%g, sum(psi_new)=%g" % (dot(conj(psi_new),dot(Ht,psi_new)), sum(abs(psi_new)**2))

    #print 
    line.set_data(p,abs(psi_new))
    time_text.set_text(time_template % (t))

    psi_old = psi_new
    
    return line, time_text

ani = animation.FuncAnimation(fig, animate, arange(1, Nt+1), interval=50, blit=True, init_func=init,repeat=False)
plt.show()
#plt.plot(linspace(0,tfinal,Nt+1),energies_t)
"""