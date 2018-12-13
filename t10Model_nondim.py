# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 17:19:53 2018

@author: john
"""
import numpy as np
import matplotlib.pyplot as plt; plt.close('all')
import scipy.sparse as sparse
import time as time
import scipy.sparse.linalg as linalg
import matplotlib.animation as animation
from finiteDifferencing import toolbox as fdt
#from matplotlib.animation import FuncAnimation as _FuncAnimation

#import hbtepLib as hbt; reload(hbt)
#import johnsUnderDevelopementToolbox as john; reload(john); 


########################################
### sub-functions

#import numpy as _np


 

def findDomainRanges(r,r_surf_index,W,oldMiddleRange):
	"""
	Because the islands grow and shrink, the location of the boundaries are in 
	constant flux.  This code figures out the boundaries and returns the 
	indices associates with all three regions.  It also returns a boolean 
	corresponding to if any of the three regions changed in size
	"""
	
	innerBCIndex=findNearest(r,r[r_surf_index]-W/2)
	if innerBCIndex==r_surf_index:
		innerBCIndex-=1;
		
	outerBCIndex=findNearest(r,r[r_surf_index]+W/2)
	if outerBCIndex==r_surf_index:
		outerBCIndex+=1;
		
	innerIndexRange=range(0,innerBCIndex+1)
	middleIndexRange=range(innerBCIndex+1,outerBCIndex)
	outerIndexRange=range(outerBCIndex,len(r))
	
	if len(oldMiddleRange)!=len(middleIndexRange):
		domainChange=True
	else:
		domainChange=False
		
	return (innerIndexRange,middleIndexRange,outerIndexRange,domainChange)
		

def calcCurrentProfileFromIP(r,r_limiter,radialFunction,params,iP,j0GuessLeft=1e5,j0GuessRight=1e7,j0Guess=1e6,errorTol=1e-6,verbose=False):
	"""
	The references only provide I_P and do not provide j(r=0).  This 
	subfunction makes a guess at j(0) and calculates j(r) with the provided
	q-profile function.  It then iterates until the integral is equal to IP.  
	
	Parameters
	----------
	r : numpy.array
		radial coordinate array
	r_limiter : float
		radial location of the limiter
	radialFunction : function(r,params)
		returns radial density current distribution
	params : list
		list of parameters to pass to radialFunction
	iP : float
		plasma current [amps]
	j0GuessLeft : float
		lower bound of j0 guess value
	j0GuessRight : float
		upper bound of j0 guess value
	errorTol : float
		error tolerance to end iteration
		
	Return
	------
	j : np.array
		radial current density where it's intergal = iP
	
	References
	----------
	http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
	"""
	
#	j0Guess=np.mean([j0GuessLeft,j0GuessRight])
	j=radialFunction(r,[j0Guess,params[1],params[2]])
	ITotal=fdt.firstOrderIntegration(r,j*r)*2*np.pi
	error=(ITotal-iP)/iP
	
	count=0
	if verbose==True:
		print('Starting iterative solver to calculated current density given the plasma current')
	while(np.abs(error)>errorTol):
		count+=1
	
		if error<0:
			j0GuessLeft=j0Guess
		else:
			j0GuessRight=j0Guess
			
		j0Guess=np.mean([j0GuessLeft,j0GuessRight])
		
		j=radialFunction(r,[j0Guess,params[1],params[2]])
		
		ITotal=fdt.firstOrderIntegration(r,j*r)*2*np.pi
		error=(ITotal-iP)/iP
		if verbose==True:
			print('count: %d, \t error: %.6f \t guess: %.3f, \t I: %.1f' % (count,error,j0Guess,ITotal))

	return j
	
def findNearest(array,value):
	"""
	search through `array` and returns the `index` of the cell closest to the 
	`value`.   `array` should be sorted in ascending order
	
	Parameters
	----------
	array : numpy.array
		data array to search through
	value : float (or int)
		value to look for in array
		
	Return
	------
	index : int
		index of value in array that is closest to value
	
	References
	----------
	http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
	"""
	index = (np.abs(array-value)).argmin()
	# value = array[index] 
	return index 
	# return index, value # uncomment to return both the index AND the value
	
	
	
## misc parameter calculations
	
def calcBeta(rho,Gamma,rho_limiter,rho_limiter_index,midRange,chiC,chiS):
#	mu0=4*np.pi*1e-7
	drho=rho[1]-rho[0]
	
	betaC=np.zeros(len(r))
	betaS=np.zeros(len(r))
	
	# add current source term
	# note that the discretized delta function requires dividing by drho
	betaC[rho_limiter_index]=-rho_limiter*Gamma/drho

	# impose boundary conditions
	betaC[0]=0
	betaC[-1]=0
	betaC[midRange[0]-1]=chiC*psi1/psi0
	betaC[midRange]=chiC*psi1/psi0
	betaC[midRange[-1]+1]=chiC*psi1/psi0
	betaS[0]=0
	betaS[-1]=0
	betaS[midRange[0]-1]=chiS*psi1/psi0
	betaS[midRange]=chiS*psi1/psi0
	betaS[midRange[-1]+1]=chiS*psi1/psi0
	
	return (betaC,betaS)

	
	
def createA(r,gamma1,gamma0,gammaM1):
	# calculate A matrix
	A=createTriDiag(gamma1,gamma0,gammaM1)
	
	A[0,0]=1							# enforces left BC
	A[0,1]=0							# enforces left BC
	A[-1,-1]=1						  # enforces right BC
	A[-1,-2]=0						  # enforces right BC
	
#	return sparse.dia_matrix(A)
	return sparse.csc_matrix(A) 
#	return sparse.csr_matrix(A) 
	
	

	
	
## finite differencing codes
	
#def firstOrderIntegration(x,y):
#	dx=x[1]-x[0]
#	return np.sum(dx*y)

#def firstOrderCenterDiff(x,y):
#	# 1st order center difference
#	dx=x[1]-x[0]
#	dydx=np.zeros(len(x))
#	dydx[0]=(y[1]-y[0])/dx
#	dydx[-1]=(y[-1]-y[-2])/dx
#	for i in range(1,len(x)-1):
#		dydx[i]=(y[i+1]-y[i-1])/(2*dx)
#		
#	return dydx
	
def createTriDiag(diag1,diag2,diag3):
	# tri-diagonal matrix
	A=np.zeros((len(diag1),len(diag1)))
	A[0,0]=diag2[0];
	A[0,1]=diag3[0];
	for i in range(1,len(diag1)-1):
		A[i,i-1]=diag1[i]
		A[i,i  ]=diag2[i]
		A[i,i+1]=diag3[i]
	A[-1,-2]=diag1[-1]
	A[-1,-1]=diag2[-1]
	return A
	
	
## current profiles
	
	
def wessonCurrentModel(r,params):
	# 
	# params = [1,0.27,3]
	j0=params[0] # j(r=0)=j0
	r0=params[1] # plasma edge (last closed flux surface)
#	q_surf=params[2] # q value at dominant surface
#	l=q_surf-1
	l=params[2]
	j=j0*(1-(r/r0)**2)**(l)
	j[np.where(r>r0)]=0
	return j
	
   

## q-profiles   
   
def quadraticQProfile(r,q0,r1,q1):
	"""
	Fit a quadratic function to the provided BCs to get q(r)
	"""
	# quadratic model, q ~ r**2
	# q(r=0) and q1(r=r1) are inputs
	c=(q1-q0)/r1**2;
	q=c*r**2+q_offset
	return q
	
def cylindricalQApproximation(r,r_limiter,l):
	"""
	Recommended in Ivanov's 2014 paper. 
	
	Notes
	-----
	The original source for q(r) is only valid for r<=a.  
	To correct for this, I solved \int B_{\theta} dl = \mu I_p and 
	q=\frac{rB_z}{RB_{\theta}} to provide q all the way out to r=b.
	"""
	q=  2*(l+1)*BT/(mu0*j[0]*R)*(r/r_limiter)**2/(1-(1-(r/r_limiter)**2)**(l+1))
	q[0]=q[1]
	i=np.where(q>0)[0]
	for k in range(i[-1]+1,len(q)):
		q[k]=2*np.pi*r[k]**2*BT/(R*mu0*iP)
	return q
	
def qProfileModel(r,j,BT,R):
	mu0=4*np.pi*1e-7
	q=np.zeros(len(r))
	for i in range(0,len(r)):
		q[i]=1/np.average(j[0:i+1])
	q*=2*BT/(mu0*R)
	return q
	
	
## plots
	
def plotInitialConditions(y2Axis=False):
	# plot j(r) 
	f,axx = plt.subplots(2,sharex=True)
	ax=axx[0]
	p1=ax.plot(rho*r0,j,'k',label='current profile')
	#ax.set_xlabel('minor radius (m)')
	ax.set_ylabel(r'current density (A/m$^2$)')
	ylim=ax.get_ylim()
	p3=ax.plot((r_surf,r_surf),ylim,'--',label=r'r$_{surf}$')
	p4=ax.plot((r_limiter,r_limiter),ylim,'--',label=r'r$_{limiter}$')
	p5=ax.plot((r_wall,r_wall),ylim,'--',label=r'r$_{wall}$')
	ax.set_ylim(ylim)
	
	
	# optional  dj(r)/dr plot  
	if y2Axis==True:
		ax2=ax.twinx()
		p2=ax2.plot(rho*r0,djdrho/r0,'r',label='current profile derivative')
		ax2.set_ylabel(r'current density derivative (A/m$^3$)',color='r')
		lns = p1+p2+p3+p4+p5
	else:
		lns=p1+p3+p4+p5
	labs = [i.get_label() for i in lns]
	ax.legend(lns, labs)#, loc=0)
	ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	
	# plot q(r) 
	ax=axx[1]
	p1=ax.plot(rho*r0,q,'k',label='q(r)')
	ax.set_xlabel('minor radius (m)')
	ax.set_ylabel(r'q')
	ylim=ax.get_ylim()
	p3=ax.plot((r_surf,r_surf),ylim,'--',label=r'r$_{surf}$')
	p4=ax.plot((r_limiter,r_limiter),ylim,'--',label=r'r$_{limiter}$')
	p5=ax.plot((r_wall,r_wall),ylim,'--',label=r'r$_{wall}$')
	ax.set_ylim(ylim)
	#ax.legend()
	if y2Axis==True: # opertional, also plot deriv of 1/q
		ax2=ax.twinx()
		p2=ax2.plot(rho*r0,dmudrho/r0,'r',label=r'$\frac{\partial (1/q)}{\partial r}$')
		ax2.set_ylabel(r'$\frac{\partial (1/q)}{\partial r}$',color='r')
		lns = p1+p2+p3+p4+p5
	else:
		lns = p1+p3+p4+p5	
	labs = [i.get_label() for i in lns]
	ax.legend(lns, labs)#, loc=0)
	
def plotFinalState(tStart=None,tStop=None,title=''):
	if tStart==None:
		tStart=t[0]
	if tStop==None:
		tStop=t[-1]
	iStart=findNearest(t,tStart)
	iStop=findNearest(t,tStop)
	f, axarr = plt.subplots(3, sharex=True)
	axarr[1].plot(t0*tau[iStart:iStop+1],BC[iStart:iStop+1]*1e4,'r',label=r'B$_C(r_{wall})$')
	axarr[1].plot(t0*tau[iStart:iStop+1],BS[iStart:iStop+1]*1e4,'b',label=r'B$_S(r_{wall})$')
	axarr[1].set_ylabel('Gauss')
	axarr[1].legend()
	axarr[2].plot(t0*tau[iStart:iStop+1],W[iStart:iStop+1]*W0,'r',label='island width')
	axarr[2].set_ylabel('m')
	axarr[2].legend()
	axarr[0].plot(t0*tau[iStart:iStop+1],Gamma[iStart:iStop+1],label='Sourced Current')
	axarr[0].set_xlabel('Time (s)')
	axarr[0].set_ylabel('A')
	axarr[0].legend()
	axarr[2].set_xlim([t0*tau[iStart],t0*tau[iStop]])
	axarr[0].set_title(title)

def psiFrame(i):
#	psi0=psi1
	fig,ax = plt.subplots()
	ax2=ax.twinx()  
	p1=ax.plot(rho,phiC[:,i]*psi0,label=r'$\psi_C$')
	p2=ax.plot(rho,phiS[:,i]*psi0,'--',label=r'$\psi_S$')
	p3=ax2.plot(rho[inRange],betaC[inRange,i]*psi0,'r',label=r'$\beta_C$')   
	ax2.plot(rho[outRange],betaC[outRange,i]*psi0,'r')   
	lns = p1+p2+p3
	labs = [count.get_label() for count in lns]
	ax.legend(lns, labs)
	ax.set_ylim([-0.0002,0.0002])
	ax.set_ylabel(r'$\psi$')
	ax2.set_ylabel(r'$\beta_C$',color='r')
	ax2.set_ylim([-0.0002,0.0002]) 
	ax.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
	ax2.ticklabel_format(style='sci',axis='y',scilimits=(0,0)) 

class animatePlot(object):
	"""An animated scatter plot using matplotlib.animations.FuncAnimation."""
	def __init__(self,step=10):
		
		
		self.fig, self.ax = plt.subplots()
		self.ax2=self.ax.twinx()
		
		# initialize
		i=0
		self.p1=self.ax.plot(rho,phiC[:,i]*psi0,label=r'$\Psi_C$') #,animated=True
		self.p2=self.ax.plot(rho,phiS[:,i]*psi0,'--',label=r'$\Psi_S$')
		self.p3=self.ax2.plot(rho,betaC[:,i]*psi0,'r',label=r'$\beta_C$') 
		lns = self.p1+self.p2+self.p3
		labs = [count.get_label() for count in lns]
		self.ax.legend(lns, labs)
		self.ax.set_ylim([-0.00015,0.00015])
		self.ax2.set_ylabel(r'$\beta_C$',color='r')
		self.ax2.set_ylim([-0.00015,0.00015]) 
		self.ax.set_title('Time = %.6f/%.6f' % (t[i], t[-1]))
		
		
		self.ani = animation.FuncAnimation(self.fig, self._update, frames=np.arange(0,len(t),step),interval=200) # , blit=True


	def _update(self, i):
		"""Update the plot."""
		
		self.p1[0].set_ydata(phiC[:,i]*psi0)
		self.p2[0].set_ydata(phiS[:,i]*psi0)
#		self.p3[0].set_xdata(r[inRange])
		self.p3[0].set_ydata(betaC[:,i]*psi0)
#		self.p3[0].set_ydata(betaC[inRange,i])
#		self.p4[0].set_xdata(r[outRange])
#		self.p4[0].set_ydata(betaC[outRange,i])
		self.ax.set_title('Time = %.6f/%.6f' % (t[i], t[-1]))
		return self.p1,

	def show(self):
		plt.show()
		
	def saveAsGif(self,fileName,dpi=75):		
		self.ani.save(fileName, dpi=dpi, writer='imagemagick')
	
		
########################################
### inputs/constants
	
## inputs
nPoints=1000 +1    	# number of radial grid points
dt=1e-5	          	# time step [seconds]

J0=200#200	          	# sourced current amplitude [Amps].  (Not for feedback)
fbGain=20          	# feedback gain absolute value  (feedback only)

# uncomment one of the following operating modes
#operatingMode="step"
#operatingMode="feedforward"
#operatingMode="feedback_suppression"
operatingMode="feedback_amplification"
#operatingMode="custom"
#operatingMode="noCurrent"

# uncomment machine 
machine='T10'
#machine='HBT'

## physical constants
mu0=4*np.pi*1e-7
		
		
########################################
### main code

## machine constants
if machine=='T10':
	m=2
	n=1
	R=1.5
	BT=2.5
	iP=250e3
	Omega=1e3*2*np.pi # mode frequency
	omegaR=1/0.01 # default 1/.01. Note that 1/.1 results in 5 gauss modes
	k=np.pi
	r_wall=.39
	r_limiter=0.27
	q_offset=0.7/.85 #q_offset and q_limiter appear to have very little to do with actual q values in the q profile....
	q_limiter=2.4
	
	psiC_s_guess=2e-4  	# guess at \Psi_C initial value at resonant surface
	psiS_s_guess=1e-4  	# guess at \Psi_S initial value at resonant surface

elif machine=='HBT':
	m=2
	n=1
	R=.92
	BT=.35
	iP=10e3
	Omega=8e3*2*np.pi # mode frequency
	omegaR=1/.001
	k=np.pi
	r_wall=.16
	r_limiter=0.15
	q_offset=.9 #q_offset and q_limiter appear to have very little to do with actual q values in the q profile....
	q_limiter=3

	
	psiC_s_guess=2e-5  	# guess at \Psi_C initial value at resonant surface
	psiS_s_guess=1e-5  	# guess at \Psi_S initial value at resonant surface
	dt=.1e-5
	

# variable/parameter non-dimensional constants (the remaining terms are defined later)
r0=r_limiter
J1=1
psi0=mu0*J1*m/2
t0=1./Omega

# create radial domain
r=np.linspace(0,r_wall,nPoints)
dr=r[1]-r[0]

# create non-dim radial domain
rho=r/r0
drho=rho[1]-rho[0]

# create time domain
if operatingMode=="step":
	tStop=30e-3
elif operatingMode=="noCurrent":
	tStop=35e-3
elif operatingMode=="feedforward":
	tStop=55e-3
elif operatingMode=="feedback_suppression" or operatingMode=="feedback_amplification":
	tStop=35e-3
elif operatingMode=="custom":
	tStop=35e-3
t=np.arange(0,tStop+dt,dt)

# create non-dim time domain
tau=t/t0
dTau=tau[1]-tau[0]

# create figure title
title=operatingMode+'. N=%d. dt=%1.1e.'%(nPoints,dt) 

# init sourced current 
J=np.zeros(len(t))

# operating mode and set currents
if operatingMode=="step":
	J[np.where(t>t[-1]/2)]=J0
	feedback=False
elif operatingMode=="feedforward":
	J0=200
	index1=np.where((t>=1e-2)&(t<=2.5e-2))[0]
	J[index1]=J0*np.sin(2*np.pi*1.5e3*(t[index1]-t[index1][0]))
	index2=np.where((t>=3e-2)&(t<=3.5e-2))[0]
	J[index2]=J0*np.sqrt(1-(t[index2]-0.5e-2-3e-2)**2/(.5e-2)**2)
	index3=np.where((t>=3.5e-2)&(t<=4e-2))[0]	
	J[index3]=J0
	index4=np.where((t>=4e-2)&(t<=4.5e-2))[0]
	J[index4]=J0*np.sqrt(1-(t[index4]-4e-2)**2/(.5e-2)**2)
	feedback=False
elif operatingMode=="feedback_suppression" or operatingMode=="feedback_amplification":
	feedback=True
	if operatingMode=="feedback_suppression":
		fbGain=np.abs(fbGain)
	else:
		fbGain=-np.abs(fbGain)
	title=title+' Gain=%.2e'%(fbGain)
	timeFeedbackOn=15e-3
	timeFeedbackOff=30e-3
elif operatingMode=="custom":
	index=np.where(t>t[-1]/2)[0]
	J[index]=J0*np.cos(2*np.pi*0.5e3*t[index]+np.pi)
	feedback=False
elif operatingMode=="noCurrent":
#	index=np.where(t>t[-1]/2)[0]
#	J[index]=J0*np.cos(2*np.pi*0.5e3*t[index]+np.pi)
	feedback=False
else:
	print('No valid operating mode provided.  Stopping code.')
	raise SystemExit
	
# create non-dim current
Gamma=J/J1

# current profile and derivative profile
l=q_limiter/q_offset-1
j=calcCurrentProfileFromIP(r,r_limiter=r_limiter,iP=iP,
						   radialFunction=wessonCurrentModel, 
						   params=[1,r_limiter,l],j0Guess=2783578.873)
#djdr=firstOrderCenterDiff(r,j)
djdrho=fdt.firstOrderSingleCenterDiff(rho,j)

# create q profile
#q=quadraticQProfile(r,q0=q_offset,r1=r_limiter,q1=q_limiter)
q=cylindricalQApproximation(r,r_limiter,l)

# calculate matrix bands (diagonal terms)
gamma1=(1./(2*drho)+rho/drho**2)#calcGamma1(rho)
gamma0= -2.*rho/drho**2-m**2/rho-mu0*R/float(BT)*djdrho/(1/q-float(n)/m)#calcGamma0(rho,djdrho,q)
if rho[0]==0:
	gamma0[0]=gamma0[1]
gammaM1=(-1./(2*drho)+rho/drho**2)#calcGammaM1(rho)

# find rational surface
rho_surf_index=findNearest(q,float(m)/n)
r_surf=r[rho_surf_index]
rho_surf=rho[rho_surf_index]

# limiter location and index
rho_limiter=r_limiter/r0
#r_limiter_index=findNearest(r,r_limiter)
rho_limiter_index=findNearest(rho,rho_limiter)

# calculate mu, its radial derivative, and its value at the mode surface
#mu=1/q
#dmudr=firstOrderCenterDiff(r,1./q)
dmudrho=fdt.firstOrderSingleCenterDiff(rho,1./q)
#dmudr_surf=dmudr[rho_surf_index]
dmudrho_surf=dmudrho[rho_surf_index]

# psi1
psi1=(-psi0**2 * t0**2 * k**2 * r_limiter**2 * omegaR**2 * rho_surf * BT * dmudrho_surf / (16*R))	**(1./3.)

# W0
W0=np.sqrt(- 16 * R * psi1 / (rho_surf * BT * dmudrho_surf))

# initialize beta
betaC=np.zeros((len(r),len(t)))
betaS=np.zeros((len(r),len(t)))

# initialize island width
W=np.zeros(len(t))

# initialize magnetic field measurements at r_wall
BC=np.zeros(len(t))
BS=np.zeros(len(t))

# initialize PsiC and PsiS
phiC=np.zeros((len(rho),len(tau)))
phiS=np.zeros((len(rho),len(tau)))

# initialize PsiC and PsiS at the surface
chiC=np.zeros(len(tau))
chiC[0]=psiC_s_guess/psi1
chiS=np.zeros(len(tau))
chiS[0]=psiS_s_guess/psi1

# set reference timer
timerRef=time.time()

# tracks how often the domains are resized
domainChange=np.zeros(len(t),dtype=bool)

# main loop
#iStop=len(t)
for i in range(0,len(t)):#len(t)):
#	print i
	
	# update non-dim island width
	W[i]=(chiC[i]**2+chiS[i]**2)**(0.25)
	
	# break up domain into inner (r<r_surface-W/2), outer (r>r_surface+W/2), and
	# middle (r_surface-W/2 <= r <= r_surface+W/2)
	if i == 0:
		midRange=[]
  	(inRange,midRange,outRange,domainChange[i])=findDomainRanges(rho,rho_surf_index,W[i]*W0/r0,midRange)
  
	# optional feedback
	if feedback==True:
		if t[i]>timeFeedbackOn and t[i]<timeFeedbackOff:
			Gamma[i]=fbGain*(BS[i-1]-BS[i-2])/dTau/J1/t0 # backward euler derivative of B_S
		
	# update RHS (beta) terms
	(betaC[:,i],betaS[:,i])=calcBeta(rho,Gamma[i],rho_limiter,rho_limiter_index,midRange,chiC[i],chiS[i])
	
	# create matrices
	if domainChange[i]:
		AInner=createA(rho[inRange],gamma1[inRange],gamma0[inRange],gammaM1[inRange])
		AOuter=createA(rho[outRange],gamma1[outRange],gamma0[outRange],gammaM1[outRange])

	# solve BVP
	phiC[inRange,i]=linalg.spsolve(AInner,betaC[inRange,i])
	phiS[inRange,i]=linalg.spsolve(AInner,betaS[inRange,i])
	phiC[outRange,i]=linalg.spsolve(AOuter,betaC[outRange,i])
	phiS[outRange,i]=linalg.spsolve(AOuter,betaS[outRange,i])
	phiC[midRange,i]=chiC[i]*psi1/psi0
	phiS[midRange,i]=chiS[i]*psi1/psi0
	
	# solve for field at r=b
	BC[i]=psi0/r0*(phiC[outRange[-1],i]-phiC[outRange[-2],i])/drho
	BS[i]=psi0/r0*(phiS[outRange[-1],i]-phiS[outRange[-2],i])/drho
	
	# solve for \Delta'
	deltaP_C_nondim=((phiC[midRange[-1]+2,i]-phiC[midRange[-1]+1,i])/drho-(phiC[midRange[0]-1,i]-phiC[midRange[0]-2,i])/drho)
	deltaP_S_nondim=((phiS[midRange[-1]+2,i]-phiS[midRange[-1]+1,i])/drho-(phiS[midRange[0]-1,i]-phiS[midRange[0]-2,i])/drho)

	# evolve in time - forward Euler
	if i < len(tau)-1:
		chiC[i+1]=chiC[i]+ dTau*( deltaP_C_nondim/W[i] - chiS[i])
		chiS[i+1]=chiS[i]+ dTau*( deltaP_S_nondim/W[i] + chiC[i])

	# print progress
	if (time.time()-timerRef)>10: # print status after every 10 seconds
		print("step=%d/%d, \t time=%.6f" % (i,len(t),t[i]))
		timerRef=time.time()
		

# plot initial conditions
#plotInitialConditions(y2Axis=True)
plotInitialConditions(y2Axis=False)

# plot final state
plotFinalState(0.15e-2,t[-1],title)
#plotFinalState()

# display the percentage of times that the A matrices needed to be recreated
temp=np.average(domainChange)*1e2
print('Domains changed %.2f %% of the time' % temp)
#plt.figure()
#plt.plot(t,domainChange)

# create animation plot
if False:
	a = animatePlot()
	a.show()
	a.saveAsGif('animation.gif')

