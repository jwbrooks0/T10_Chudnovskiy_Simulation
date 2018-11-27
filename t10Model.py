# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 17:19:53 2018

@author: john
"""
import numpy as np
import matplotlib.pyplot as plt; plt.close('all')
import scipy.sparse as sparse
import time as time
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
    indices associates with all three regions.  
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
        

def calcCurrentProfileFromIP(r,r_limiter,radialFunction,params,iP,j0GuessLeft=1e5,j0GuessRight=1e7,errorTol=1e-6):
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
    
    j0Guess=np.mean([j0GuessLeft,j0GuessRight])
    j=radialFunction(r,[j0Guess,params[1],params[2]])
    ITotal=firstOrderIntegration(r,j*r)*2*np.pi
    error=(ITotal-iP)/iP
    
    count=0
    
    while(np.abs(error)>errorTol):
        count+=1
    
        if error<0:
            j0GuessLeft=j0Guess
        else:
            j0GuessRight=j0Guess
            
        j0Guess=np.mean([j0GuessLeft,j0GuessRight])
        
        j=radialFunction(r,[j0Guess,params[1],params[2]])
        
        ITotal=firstOrderIntegration(r,j*r)*2*np.pi
        error=(ITotal-iP)/iP
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
    
def calcBeta(r,Jr,r_limiter,r_limiter_index,midRange,psiC_s,psiS_s):
    mu0=4*np.pi*1e-7
    dr=r[1]-r[0]
    
    betaC=np.zeros(len(r))
    betaS=np.zeros(len(r))
    
    # add current source term
    iota=m*Jr/2/r_limiter/dr
    betaC[r_limiter_index]=iota*mu0*r_limiter

    # impose boundary conditions
    betaC[0]=0
    betaC[-1]=0
    betaC[midRange[0]-1]=psiC_s
    betaC[midRange]=psiC_s
    betaC[midRange[-1]+1]=psiC_s
    betaS[0]=0
    betaS[-1]=0
    betaS[midRange[0]-1]=psiS_s
    betaS[midRange]=psiS_s
    betaS[midRange[-1]+1]=psiS_s
    
    return (betaC,betaS)

    
def calcAlpha(r,djdr,q):
    return m**2/r+mu0*R/float(BT)*djdr/(1/q-float(n)/m)

def calcGamma1(r):
    dr=r[1]-r[0]
    return (1/(2*dr)+r/dr**2)
def calcGamma0(r,djdr,q):
    alpha=calcAlpha(r,djdr,q)
    dr=r[1]-r[0]
    return -2*r/dr**2-alpha
def calcGammaM1(r):
    dr=r[1]-r[0]
    return (-1/(2*dr)+r/dr**2)
    
def width(r_surf,dmudr_surf,psiC_s,psiS_s):
    return 4*np.sqrt(np.sqrt(psiC_s**2+psiS_s**2)/(-r_surf*BT*dmudr_surf/R))

def calcDeltaPrime(dr,psiA,psiB):
    dPsiA=(psiA[-1]-psiA[-2])/dr 
    dPsiB=(psiB[1]-psiB[0])/dr 
    return (dPsiB-dPsiA)/psiB[0]

def createA(r,gamma1,gamma0,gammaM1):
    # calculate A matrix
    A=createTriDiag(gamma1,gamma0,gammaM1)
    
    A[0,0]=1                            # enforces left BC
    A[0,1]=0                            # enforces left BC
    A[-1,-1]=1                          # enforces right BC
    A[-1,-2]=0                          # enforces right BC
    
    return sparse.dia_matrix(A)
    
    

    
    
## finite differencing codes
    
def firstOrderIntegration(x,y):
    dx=x[1]-x[0]
    return np.sum(dx*y)

def firstOrderCenterDiff(x,y):
    # 1st order center difference
    dx=x[1]-x[0]
    dydx=np.zeros(len(x))
    dydx[0]=(y[1]-y[0])/dx
    dydx[-1]=(y[-1]-y[-2])/dx
    for i in range(1,len(x)-1):
        dydx[i]=(y[i+1]-y[i-1])/(2*dx)
        
    return dydx
    
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
#    q_surf=params[2] # q value at dominant surface
#    l=q_surf-1
    l=params[2]
    j=j0*(1-(r/r0)**2)**(l)
    j[np.where(r>r0)]=0
    return j
    
   

## q-profiles   
   
def quadraticQProfile(r,q0,r1,q1):
    # quadratic model, q ~ r**2
    # q(r=0) and q1(r=r1) are inputs
    c=(q1-q0)/r1**2;
    q=c*r**2+q_offset
    return q
    
def cylindricalQApproximation(r,r_limiter,l):
    """
    Recommended in Ivanov's 2014 paper. However, it is only valid for r<=a.  
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
    
    
#def addBoundariesToPlot(ax,r_surf,r_limiter,r_wall):
#    ylim=ax.get_ylim()
#    ax.plot((r_surf,r_surf),ylim,'--',color='grey')
#    ax.plot((r_limiter,r_limiter),ylim,'--',color='grey')
#    ax.plot((r_wall,r_wall),ylim,'--',color='grey')
#    ax.set_ylim(ylim)
    
def plotInitialConditions():
    
    
    # plot j(r) and dj(r)/dr    
    f,axx = plt.subplots(2,sharex=True)
    ax=axx[0]
    ax2=ax.twinx()
    p1=ax.plot(r,j,'k',label='current profile')
    p2=ax2.plot(r,djdr,'r',label='current profile derivative')
    #ax.set_xlabel('minor radius (m)')
    ax.set_ylabel(r'current density (A/m$^2$)')
    ax2.set_ylabel(r'current density derivative (A/m$^3$)',color='r')
    ylim=ax.get_ylim()
    p3=ax.plot((r_surf,r_surf),ylim,'--',label=r'r$_{surf}$')
    p4=ax.plot((r_limiter,r_limiter),ylim,'--',label=r'r$_{limiter}$')
    p5=ax.plot((r_wall,r_wall),ylim,'--',label=r'r$_{wall}$')
    ax.set_ylim(ylim)
    #ax.legend()
    lns = p1+p2+p3+p4+p5
    labs = [i.get_label() for i in lns]
    ax.legend(lns, labs)#, loc=0)
    
    
    # plot q(r) and d(1/q(r))/dr
    #f,ax = plt.subplots(1)
    ax=axx[1]
    ax2=ax.twinx()
    p1=ax.plot(r,q,'k',label='q(r)')
    p2=ax2.plot(r,dmudr,'r',label=r'$\frac{\partial (1/q)}{\partial r}$')
    ax.set_xlabel('minor radius (m)')
    ax.set_ylabel(r'q')
    ax2.set_ylabel(r'$\frac{\partial (1/q)}{\partial r}$',color='r')
    ylim=ax.get_ylim()
    p3=ax.plot((r_surf,r_surf),ylim,'--',label=r'r$_{surf}$')
    p4=ax.plot((r_limiter,r_limiter),ylim,'--',label=r'r$_{limiter}$')
    p5=ax.plot((r_wall,r_wall),ylim,'--',label=r'r$_{wall}$')
    ax.set_ylim(ylim)
    #ax.legend()
    lns = p1+p2+p3+p4+p5
    labs = [i.get_label() for i in lns]
    ax.legend(lns, labs)#, loc=0)
    
def plotFinalState(tStart,tStop):
    iStart=findNearest(t,tStart)
    iStop=findNearest(t,tStop)
    f, axarr = plt.subplots(3, sharex=True)
    axarr[0].plot(t[iStart:iStop+1],BC[iStart:iStop+1]*1e4,'r',label=r'B$_C(r_{wall})$')
    axarr[0].plot(t[iStart:iStop+1],BS[iStart:iStop+1]*1e4,'b',label=r'B$_S(r_{wall})$')
    axarr[0].set_ylabel('Gauss')
    axarr[0].legend()
    axarr[1].plot(t[iStart:iStop+1],W[iStart:iStop+1],'r',label='island width')
    axarr[1].set_ylabel('m')
    axarr[1].legend()
    axarr[2].plot(t[iStart:iStop+1],J[iStart:iStop+1],label='Sourced Current')
    axarr[2].set_xlabel('Time (s)')
    axarr[2].set_ylabel('A')
    axarr[2].legend()
    axarr[2].set_xlim([0,tStop])

def psiFrame(i):
    fig,ax = plt.subplots()
    ax2=ax.twinx()  
    p1=ax.plot(r,PsiC[:,i],label='PsiC')
    p2=ax.plot(r,PsiS[:,i],'--',label='PsiS')
    p3=ax2.plot(r[inRange],betaC[inRange,i],'r',label=r'$\beta_C$')   
    ax2.plot(r[outRange],betaC[outRange,i],'r')   
    lns = p1+p2+p3
    labs = [count.get_label() for count in lns]
    ax.legend(lns, labs)
    ax.set_ylim([-0.0002,0.0002])
    ax2.set_ylabel(r'$\beta_C$',color='r')
#    ylim=ax2.get_ylim()
    ax2.set_ylim([-0.0002,0.0002])  
        
        
########################################
### inputs/constants
    
## inputs
nPoints=1000 +1       # number of radial grid points
tStop=5.5e-2#1e-2 #5.5e-2         # duration of simulation [seconds]
dt=1e-5           # time step [seconds]

J0=200#200            # sourced current [Amps]
phi0=0*np.pi

psiC_s_guess=2e-4#e-7         # guess at \Psi_C initial value at resonant surface
psiS_s_guess=psiC_s_guess/2.#1e-6#e-7         # guess at \Psi_S initial value at resonant surface

## machine constants
m=2
n=1
R=1.5
BT=2.5
iP=250e3
Omega=1e3*2*np.pi
omegaR=1/0.01
k=np.pi
r_wall=.39
r_limiter=0.27
q_offset=0.7/.85
q_limiter=2.4

## physical constants
mu0=4*np.pi*1e-7
        
        
########################################
### main code
        
# create radial domain
r=np.linspace(0,r_wall,nPoints)
dr=r[1]-r[0]

# create time domain
t=np.arange(0,tStop+dt,dt)

# derived constants
zeta=dt*k*r_limiter**2*omegaR
eta=dt*Omega

# create sourced current 
J=np.zeros(len(t))#np.sin(Omega*t+phi0)
if False:
    J[np.where(t>t[-1]/2)]=J0
if True:
    index1=np.where((t>=1e-2)&(t<=2.5e-2))[0]
    J[index1]=J0*np.sin(2*np.pi*1.5e3*(t[index1]-t[index1][0]))
    index2=np.where((t>=3e-2)&(t<=3.5e-2))[0]
    J[index2]=J0*np.sqrt(1-(t[index2]-0.5e-2-3e-2)**2/(.5e-2)**2)
    index3=np.where((t>=3.5e-2)&(t<=4e-2))[0]    
    J[index3]=J0
    index4=np.where((t>=4e-2)&(t<=4.5e-2))[0]
    J[index4]=J0*np.sqrt(1-(t[index4]-4e-2)**2/(.5e-2)**2)

# current profile and derivative profile
l=q_limiter/q_offset-1
j=calcCurrentProfileFromIP(r,r_limiter=r_limiter,iP=iP,
                           radialFunction=wessonCurrentModel, 
                           params=[1,r_limiter,l])
djdr=firstOrderCenterDiff(r,j)

# create q profile
#q=quadraticQProfile(r,q0=q_offset,r1=r_limiter,q1=q_limiter)
q=cylindricalQApproximation(r,r_limiter,l)

# calculate gamma terms
gamma1=calcGamma1(r)
gamma0=calcGamma0(r,djdr,q)
gammaM1=calcGammaM1(r)

# find rational surface
r_surf_index=findNearest(q,float(m)/float(n))
r_surf=r[r_surf_index]

# find limiter 
r_limiter_index=findNearest(r,r_limiter)

# calculate mu, its radial derivative, and its value at the mode surface
#mu=1/q
dmudr=firstOrderCenterDiff(r,1./q)
dmudr_surf=dmudr[r_surf_index]


# initialize beta
betaC=np.zeros((len(r),len(t)))
betaS=np.zeros((len(r),len(t)))

# initialize island width
W=np.zeros(len(t))

# initialize magnetic field measurements at r_wall
BC=np.zeros(len(t))
BS=np.zeros(len(t))

# initialize PsiC and PsiS
PsiC=np.zeros((len(r),len(t)))
PsiS=np.zeros((len(r),len(t)))

# initialize PsiC and PsiS at the surface
psiC_s=np.zeros(len(t))
psiC_s[0]=psiC_s_guess
psiS_s=np.zeros(len(t))
psiS_s[0]=psiS_s_guess

# set reference timer
timerRef=time.time()

# tracks how often the domains are resized
domainChange=np.zeros(len(t),dtype=bool)

# main loop
for i in range(0,len(t)):#len(t)):
    
    # update island width
    W[i]=width(r_surf,dmudr_surf,psiC_s[i],psiS_s[i])
    
    # break up domain into inner (r<r_surface-W/2), outer (r>r_surface+W/2), and
    # middle (r_surface-W/2 <= r <= r_surface+W/2)
    if i == 0:
        midRange=[]
    (inRange,midRange,outRange,domainChange[i])=findDomainRanges(r,r_surf_index,W[i],midRange)
  
    # update betas
    (betaC[:,i],betaS[:,i])=calcBeta(r,J[i],r_limiter,r_limiter_index,midRange,psiC_s[i],psiS_s[i])
    
    # create matrices
    if domainChange[i]:
        AInner=createA(r[inRange],gamma1[inRange],gamma0[inRange],gammaM1[inRange])
        AOuter=createA(r[outRange],gamma1[outRange],gamma0[outRange],gammaM1[outRange])

    # solve BVP
    PsiC[inRange,i]=sparse.linalg.spsolve(AInner,betaC[inRange,i])
    PsiS[inRange,i]=sparse.linalg.spsolve(AInner,betaS[inRange,i])
    PsiC[outRange,i]=sparse.linalg.spsolve(AOuter,betaC[outRange,i])
    PsiS[outRange,i]=sparse.linalg.spsolve(AOuter,betaS[outRange,i])
    PsiC[midRange,i]=psiC_s[i]
    PsiS[midRange,i]=psiS_s[i]
    
    # solve for field at r=b
    BC[i]=(PsiC[outRange[-1],i]-PsiC[outRange[-2],i])/dr
    BS[i]=(PsiS[outRange[-1],i]-PsiS[outRange[-2],i])/dr
    
    # solve for \Delta'
    deltaP_C=(-(PsiC[midRange[0]-1,i]-PsiC[midRange[0]-2,i])/dr+(PsiC[midRange[-1]+2,i]-PsiC[midRange[-1]+1,i])/dr)/psiC_s[i]
    deltaP_S=(-(PsiS[midRange[0]-1,i]-PsiS[midRange[0]-2,i])/dr+(PsiS[midRange[-1]+2,i]-PsiS[midRange[-1]+1,i])/dr)/psiS_s[i]
    
    # evolve in time - forward Euler
    if i < len(t)-1:
        psiC_s[i+1]=psiC_s[i]*(1+zeta*deltaP_C/W[i])-eta*psiS_s[i]
        psiS_s[i+1]=psiS_s[i]*(1+zeta*deltaP_S/W[i])+eta*psiC_s[i]
#    psiC_s=psiC_s+dt*(k*r_limiter**2*omegaR*deltaP_C/W[i]*psiC_s-Omega*psiS_s)
#    psiS_s=psiS_s+dt*(k*r_limiter**2*omegaR*deltaP_S/W[i]*psiS_s+Omega*psiC_s)
#    psiC_s=psiC_s+(dt*k*r_limiter**2*omegaR*deltaP_C/W[i]*psiC_s-dt*Omega*psiS_s)
#    psiS_s=psiS_s+dt*(k*r_limiter**2*omegaR*deltaP_S/W[i]*psiS_s+Omega*psiC_s)
    
    # print progress
    if (time.time()-timerRef)>10:
        print("step=%d/%d, \t time=%.6f" % (i,len(t),t[i]))
        timerRef=time.time()
        
    # plot PsiC and PsiS
    if np.mod(i,400)==0:
        psiFrame(i)


# plot initial conditions
plotInitialConditions()

# plot final state
plotFinalState(0.15e-2,t[-1])
#
#

#
#plotPsiOfR(rA,rB,PsiC_A,PsiS_A,PsiC_B,PsiS_B)
#plotBeta(rA,rB,betaC_A,betaS_A,betaC_B,betaS_B)

temp=np.average(domainChange)*1e2
print('Domains changed %.2f %% of the time' % temp)
#plt.figure()
#plt.plot(t,domainChange)





#import matplotlib.pyplot as plt
import matplotlib.animation as animation
#import numpy as np

class animatePlot(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self,step=10):
        
        
        self.fig, self.ax = plt.subplots()
        self.ax2=self.ax.twinx()
        
        # initialize
        i=0
        self.p1=self.ax.plot(r,PsiC[:,i],label='PsiC') #,animated=True
        self.p2=self.ax.plot(r,PsiS[:,i],'--',label='PsiS')
        self.p3=self.ax2.plot(r,betaC[:,i],'r',label=r'$\beta_C$') 
#        self.p3=self.ax2.plot(r[inRange],betaC[inRange,i],'r',label=r'$\beta_C$') 
#        self.p4=self.ax2.plot(r[outRange],betaC[outRange,i],'r') 
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
        
        self.p1[0].set_ydata(PsiC[:,i])
        self.p2[0].set_ydata(PsiS[:,i])
#        self.p3[0].set_xdata(r[inRange])
        self.p3[0].set_ydata(betaC[:,i])
#        self.p3[0].set_ydata(betaC[inRange,i])
#        self.p4[0].set_xdata(r[outRange])
#        self.p4[0].set_ydata(betaC[outRange,i])
        self.ax.set_title('Time = %.6f/%.6f' % (t[i], t[-1]))
        return self.p1,

    def show(self):
        plt.show()
        
    def saveAsGif(self,fileName,dpi=75):        
        self.ani.save(fileName, dpi=dpi, writer='imagemagick')

#if __name__ == '__main__':
a = animatePlot()
a.show()
a.saveAsGif('temp.gif')
