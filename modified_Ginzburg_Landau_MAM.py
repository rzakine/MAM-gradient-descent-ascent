from numpy import loadtxt
from matplotlib import cm
import math
import cmath
import time
import os

import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rc

from itertools import count
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)



import time
start_time = time.time()

PI = math.pi

def normL2(h, array): #norm L2, spacing h for integration
	norm = math.sqrt( h * np.sum( np.square( array.real ) ) )
	return norm

def normL2Big(h, array): #norm L2, spacing h for integration
	norm = np.sqrt( h * np.sum( np.square( array.real ), axis=1 ) )
	return norm




##############################
""" DEFINE FUNCTIONS """
###############################

from scipy.linalg import solve_triangular
#x = solve_triangular(a, b, lower=True)

def Hamiltonian(h, rho, theta): #return an array of size Ncopy, axis=1 sums the colums
	H = np.sum( (D*(np.roll(rho,-1, axis=1) +  np.roll(rho,1, axis=1) -2* rho)/h**2 + rho -rho**3 + kappa*np.outer(np.mean(rho**2,axis=1), np.ones(L)))*theta + 0.5*aa*theta**2     ,  axis=1)
	return H

def Lagrangian(h, ds, rho, theta): #return an array of size Ncopy, axis=1 sums the colums
	#rhoDot = (np.roll(rho,-1,axis=0) - rho)/ds
	#L = h*np.sum(np.power(theta.real,2),axis=1)
	rhoDot      = (np.roll(rho,-1,axis=0) - np.roll(rho,1,axis=0))/(2*ds)
	rhoDot[0,:] = (np.roll(rho,-1,axis=0) - rho)[0,:] /(ds)
	rhoDot[Ncopy-1,:] = (-np.roll(rho,1,axis=0) + rho)[Ncopy-1,:]/(ds)
	Ham = Hamiltonian(h, rho, theta)
	L =  np.sum( rhoDot * theta, axis=1 ) - Ham
	return L






################################
""" Start MAM """
##############################

Lx = 1.
L  = 2
h  = 1.0/L
Ncopy = 400 
Tmax = 10.


s = np.linspace(0,Tmax,Ncopy)
ds  = s[1]-s[0]
dnu = s[1]-s[0]
aa = 2.   #noise amplitude


upward = True # choose if path from -1 to +1 (upward), or the opposite

D     = 0.05
kappa = 0.26

solRho = np.sort(np.roots([-1, kappa, 1, 0])) # rho-, rhos, rho+ 

dtau = 0.02

iterations = 20000
plotStep   = 2000

r = dtau/dnu


############# FIXED ARRAYS ##########

A_solve_upper_adapted = np.zeros((L,Ncopy,Ncopy))
B_solve_lower_adapted = np.zeros((L,Ncopy,Ncopy))
## Fisrt order scheme
for kk in range(0,L):
	for ii in range(0,Ncopy-1):
		A_solve_upper_adapted[kk][ii][ii]   = 1+r
		A_solve_upper_adapted[kk][ii][ii+1] = -r
	A_solve_upper_adapted[kk][Ncopy-1,Ncopy-1] = 1+dtau

for kk in range(0,L):
	for ii in range(1,Ncopy):
		B_solve_lower_adapted[kk][ii][ii]   = 1+r
		B_solve_lower_adapted[kk][ii][ii-1] = -r
	B_solve_lower_adapted[kk,0,0] = 1+dtau	


## Second order scheme
for kk in range(0,L):
	np.fill_diagonal( A_solve_upper_adapted[kk], 1. + 3*r/2. )
	# Fill up upper diagonal:
	np.fill_diagonal( A_solve_upper_adapted[kk][:-1, 1:], -2*r ) 
	# fill up 2nd upper
	np.fill_diagonal( A_solve_upper_adapted[kk][:-2, 2:], r/2. ) 
	# correct border effects
	A_solve_upper_adapted[kk][Ncopy-2,Ncopy-2] = 1 + r
	A_solve_upper_adapted[kk][Ncopy-2,Ncopy-1] = -r
	A_solve_upper_adapted[kk][Ncopy-1,Ncopy-1] = 1


	np.fill_diagonal( B_solve_lower_adapted[kk,1:Ncopy,1:Ncopy],  1. + 3*r/2. )
	# fill up lower diagonal :
	np.fill_diagonal( B_solve_lower_adapted[kk][1:, :-1], -2*r ) 
	# fill up 2nd lower
	np.fill_diagonal( B_solve_lower_adapted[kk][2:, :-2], r/2. ) 
	# correct border effects
	B_solve_lower_adapted[kk][1,1] = 1 + r
	B_solve_lower_adapted[kk][1,0] = -r
	B_solve_lower_adapted[kk,0,0] = 1




########## ARRAY CREATION

rho   = np.zeros((Ncopy,L), dtype=complex)

theta = np.zeros((Ncopy,L), dtype=complex)

U = np.zeros((Ncopy,L), dtype=complex)
U_Fourier = np.zeros((Ncopy,L), dtype=complex)

V = np.zeros((Ncopy,L), dtype=complex)
V_Fourier = np.zeros((Ncopy,L), dtype=complex)

U2_Fourier = np.zeros((Ncopy,L), dtype=complex)
V2_Fourier = np.zeros((Ncopy,L), dtype=complex)

reaction_U = np.zeros((Ncopy,L), dtype=complex)
reaction_V = np.zeros((Ncopy,L), dtype=complex)




if(upward==True):
	rho1 = solRho[0] * np.ones(L, dtype=complex)
	rho1k = np.fft.fft( rho1)
	rho2 = solRho[2] * np.ones(L, dtype=complex)
	rho2k = np.fft.fft( rho2)
else:
	rho1 = solRho[2] * np.ones(L, dtype=complex)
	rho1k = np.fft.fft( rho1)
	rho2 = solRho[0] * np.ones(L, dtype=complex)
	rho2k = np.fft.fft( rho2)


reaction_U_Fourier = np.zeros((Ncopy,L), dtype=complex)
reaction_V_Fourier = np.zeros((Ncopy,L), dtype=complex)


######  INITIAL CONDITIONS 


rho[0]       = rho1
rho[Ncopy-1] = rho2
amp = 0.4

for j in range(1,Ncopy-1):
	tt = float(j)/Ncopy
	rho[j] = rho1*(1-tt) + tt*rho2 + amp*np.square(np.sin(PI*np.arange(0,L)/L))*np.power(np.sin(PI*tt),2)




U = rho + theta
V = rho - theta



###### Evolve complex array ######
for i in range(0, iterations+1):
	

	##################################################
	'''------------ FIRST STEP: update U ---------'''
	##################################################	
	#print 'UPDATE U...'
	#Compute nonlinear term in real space
	# dtau_rho = dt_theta + dH_drho
	# dtau_theta = dt_rho - dH_dtheta
	
	dH_drho = D*(np.roll(theta,-1, axis=1) +  np.roll(theta,1, axis=1) -2* theta)/h**2 + theta - 3*rho**2 *theta + kappa*2*rho*np.outer(np.mean(theta, axis=1), np.ones(L))
	
	dH_dtheta = D*(np.roll(rho,-1, axis=1) +  np.roll(rho,1, axis=1) -2* rho)/h**2 + rho -rho**3 + kappa*np.outer(np.mean(rho**2,axis=1), np.ones(L)) + aa*theta
	
	
	reaction_U = dH_drho - dH_dtheta
	
	#Compute Fourier transform of nonlinear terms
	reaction_U_Fourier[:] = np.fft.fft( reaction_U[:] )  #size Ncopy,L

	U_Fourier[:] = np.fft.fft( U[:]) # Ncopy * L
	V_Fourier[:] = np.fft.fft( V[:])

	
	#COMPUTE EVOLUTION for mode k, temporary array (Ncopy,L)
	
	temporaryArray = U_Fourier[:,:]

	temporaryArray[Ncopy-1] = (-V_Fourier[Ncopy-1,:] + 2*rho2k[:])#implement boundary conditions

	reaction_U_Fourier[Ncopy-1,:] = 0
	
	for kk in range(0,L):
		U2_Fourier[:,kk] = solve_triangular( A_solve_upper_adapted[kk], temporaryArray[:,kk] + dtau*(reaction_U_Fourier[:,kk])) #size 2Ncopy
		


	#Compute inverse Fourier transform of UV2
	U[:] = np.fft.ifft( U2_Fourier[:]).real #size 2Ncopy,L

	rho   = 0.5*(U+V)   #size Ncopy,L
	theta = 0.5*(U-V) #size Ncopy,L

	rho[0,:] = rho1
	rho[Ncopy-1,:] = rho2
	
	#print rho[0].real
	#print rho[Ncopy-1,:].real
	U = rho + theta #size Ncopy,L
	V = rho - theta #size Ncopy,L
	
	


	##################################################
	'''------------ SECOND STEP: update V ---------'''
	##################################################	
	#print 'UPDATE V...'
	#Compute nonlinear term in real space

	
	dH_drho = D*(np.roll(theta,-1, axis=1) +  np.roll(theta,1, axis=1) -2* theta)/h**2  + (1-3*rho**2)*theta + kappa*2*rho*np.outer(np.mean(theta,axis=1), np.ones(L))
	
	dH_dtheta = D*(np.roll(rho,-1, axis=1) +  np.roll(rho,1, axis=1) -2* rho)/h**2  + rho -rho**3 + kappa*np.outer(np.mean(rho**2,axis=1), np.ones(L))  + aa*theta
	
	
	reaction_V = dH_drho + dH_dtheta


	#Compute Fourier transform of nonlinear terms
	reaction_V_Fourier[:] = np.fft.fft( reaction_V[:])
	
	U_Fourier[:] = np.fft.fft( U[:])
	V_Fourier[:] = np.fft.fft( V[:])


	
	temporaryArray = V_Fourier[:,:] 
	#implement boundary conditions
	temporaryArray[0] = (-U_Fourier[0,:] + 2*rho1k[:])
	reaction_V_Fourier[0,:] = 0
	
	
	#COMPUTE EVOLUTION for mode k
	
	for kk in range(0,L):
		V2_Fourier[:,kk] = solve_triangular(B_solve_lower_adapted[kk], temporaryArray[:,kk] + dtau*(reaction_V_Fourier[:,kk]), lower=True)  #size 2Ncopy
		

	#Compute inverse Fourier transform of UV2
	V[:] = np.fft.ifft( V2_Fourier[:]).real #size 2Ncopy,L

	rho   = 0.5*(U+V)   #size Ncopy,L
	theta = 0.5*(U-V) #size Ncopy,L


	rho[0,:]       = rho1
	rho[Ncopy-1,:] = rho2

	
	U = rho + theta #size Ncopy,L
	V = rho - theta #size Ncopy,L
	



	##################################################
	'''------------ THIRD: PLOT DATA ---------'''
	##################################################	
	    
	if( i%plotStep == 0):
		
		plt.gcf()
		fig = plt.figure(figsize=(16,5))
		fig.tight_layout(pad=4.0)
		ax  = fig.add_subplot(131)#fig.add_axes([0.1, 0.1, 0.6, 0.75])
		Lag = Lagrangian(h, dnu, rho, theta).real
		actionS = dnu* np.sum(Lag)
		
		fig.suptitle(r'$N_\mathrm{copy}=$'+str(int(Ncopy))+r', $L=$'+str(L)+r', $\Delta \tau=$'+str("%.1e"%dtau)+r', $D=$'+str("%.1e"%D) +r', $T_\mathrm{max}=$'+str("%.1f"%Tmax)+'\n'+r'$S=$'+str("%.6f"%(actionS))+', Time '+ str(i*dtau), fontsize=14)
		
		### VECTOR FIELD
		Saction = dnu* np.sum(Lagrangian(h, dnu, rho, theta).real)
		ax.scatter(rho[:,0].real, rho[:,1].real, color='darkblue', s=2.2, zorder=15)
		x = np.arange(solRho[0], solRho[2], 0.02)
		y = np.arange(solRho[0], solRho[2], 0.02)

		X, Y = np.meshgrid(x, y)
		u = D*(2*Y-2*X)/h**2  + (X - X*X*X) + kappa*(X*X+Y*Y)/2.
		v = D*(2*X-2*Y)/h**2  + (Y - Y*Y*Y) + kappa*(X*X+Y*Y)/2.

		ax.streamplot(x, y, u, v, density=1, color='grey')
		ax.set_aspect('equal')
		ax.set_xlim(solRho[0],solRho[2])
		ax.set_ylim(solRho[0],solRho[2])
		ax.set_xlabel(r'$x_1$')
		ax.set_ylabel(r'$x_2$')
		Lag = Lagrangian(h, dnu, rho, theta).real
		actionS = dnu* np.sum(Lag)

		### PLOT theta
		ax2 = fig.add_subplot(132)
	
		ax2.plot((np.arange(0,Ncopy)/(float(Ncopy-1))), theta[:,0].real, linewidth=1, label=r'$\theta_1$' )
		ax2.plot((np.arange(0,Ncopy)/(float(Ncopy-1))), theta[:,1].real, linewidth=1, label=r'$\theta_2$' )
		ax2.legend(loc='best')
		ax2.set_xlabel(r'$t/T_\mathrm{Max}$', fontsize=15)
		#ax2.set_ylabel(r'$\theta$', fontsize=15)
		ax2.set_xlim(0,1)
		ax2.set_ylim(-1,1)

		### PLOT Lagrangian and Hamiltonian
		ax3 = fig.add_subplot(133)
		ax3.plot(np.linspace(0,1,Ncopy), Lag, label=r'$L(\rho,\dot\rho)$', color='black'  ) 
		ax3.plot(np.linspace(0,1,Ncopy), Hamiltonian(h,rho, theta).real, label=r'$H(\rho,\theta)$', color='brown' ) 
		plt.legend(loc='best', fontsize=14)
		ax3.set_xlabel(r'$t/T_\mathrm{Max}$', fontsize=18)
		ax3.set_xlim(0,1)
		
		if(upward==True):
			plt.savefig('upward_L'+str(int(L))+'_N'+str(int(Ncopy))+'_D'+str(D)+'_kappa'+str(kappa)+'_dtau'+str("%.1e"%dtau)+'_'+"%09d" % (i,)+'.png', format='png', bbox_inches='tight')
		else:
			plt.savefig('downward_L'+str(int(L))+'_N'+str(int(Ncopy))+'_D'+str(D)+'_kappa'+str(kappa)+'_dtau'+str("%.1e"%dtau)+'_'+"%09d" % (i,)+'.png', format='png', bbox_inches='tight')
		plt.clf()
		plt.close()
		












