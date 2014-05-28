import math
import numpy
from scipy.special import jv
from dolfin import *
from scipy.interpolate import *

class flow_waveform():
	def __init__(self,bc,face,FC):
		#Input parameters
		self.problem=bc.problem
		self.nu=self.problem.nu
		self.R=self.problem.radius[face]
		self.T0=self.problem.T
		self.dt=self.problem.dt
		
		#Load the Fourier Coefficients
		infile_FC=open("./data/"+FC)
		self.an=[];self.bn=[]
		for line in infile_FC:
			self.an.append(float(line.split()[0]))
			self.bn.append(float(line.split()[1]))
			
		self.omega=((2.*DOLFIN_PI)/self.T0) #omega
		self.Q_mean=(389./60.) #Mean flow rate
		self.womn=self.R*(self.omega/self.nu)**0.5 #womersley number

		#Poieuille term constant
		self.poiseuille_term_constant=((2.*self.an[0]*self.Q_mean)/(numpy.pi*self.R**2.)) 

		#Pulsatile term constant
		self.pulsatile_term_numerator1=[]
		self.pulsatile_term_denominator=[]
		for n in range(1,len(self.an)):
			alpha=(n)**0.5*self.womn
			self.pulsatile_term_numerator1.append( alpha*1j**1.5*jv(0,alpha*1j**1.5) )
			self.pulsatile_term_denominator.append( alpha*1j**1.5*jv(0,alpha*1j**1.5)- 2*jv(1,alpha*1j**1.5) )

		self.r=[] #store r values
		self.pulsatile_for_all_r=[] #store the varying component of r

		#Velocity term constant
		self.velocity_term_constant=1./(numpy.pi*self.R**2.)*(self.Q_mean)
		
                #flow rate term
                self.t=[] #the flow rate series will remain the same for all 'r' values at a timestep
                self.Qn_at_current_ts=0 #store the series for current timestep for multiple r values 

		
        def flowrate_term(self,t):
                #flow rate term
                if t in self.t:
                        return self.Qn_at_current_ts #Q for all 'r' will be the same at the current timestep
                else:
                        Qn=[]
                        for n in range (1,len(self.an)):
                                Qn.append( (self.an[n]-self.bn[n]*1j)*numpy.exp(1j*n*self.omega*t))
                        self.Qn_at_current_ts=Qn #if the flow rate at current time has not been calculated then calculate it 
                        self.t.append(t) #append the correct time
                        return Qn


	def poiseuille_term(self,r):
		#Determine the Poiseuille term ofof the velocity profile
		poiseuille=self.poiseuille_term_constant* (1.-(r/self.R)**2.)
		return poiseuille 


	def pulsatile_term(self,r):
		#store the value of r
		if r in self.r:
			return self.pulsatile_for_all_r[self.r.index(r)]
		else:
			#Calculate the pulsatile term
			pulsatile=[]
			y=r/self.R
			for n in range (1,len(self.an)):
				alpha=(n)**0.5*self.womn
				Numerator=self.pulsatile_term_numerator1[n-1] - alpha*1j**1.5*jv(0,alpha*1j**1.5*y)
				pulsatile.append(Numerator/self.pulsatile_term_denominator[n-1])
			self.pulsatile_for_all_r.append(pulsatile)
			self.r.append(r)
			return pulsatile


	def velocity_term(self,r,t):
		Qn=self.flowrate_term(t)
		pulsatile=self.pulsatile_term(r)
		velocity=0
		for i in range (0,len(self.an)-1):
			velocity=velocity+self.velocity_term_constant*Qn[i]*pulsatile[i]
		velocity=numpy.real(velocity)+self.poiseuille_term(r)
		return velocity
			
        def validate_Q(self,t):
                Qn=self.flowrate_term(t)
                Q_temp=self.an[0]*self.Q_mean
                for j in range (0,len(Qn)):
                        Q_temp=Q_temp+Qn[j]*self.Q_mean
                return Q_temp.real
	


class interpolate_wom():
	def __init__(self,bc,face,FC):
		#call the womersley function
		t=bc.problem.t
		T=bc.problem.T
		R=bc.problem.radius[face]
		wom_func=flow_waveform(bc,face,FC)
		
		#Create a 2D spline of Womersley
		xres=50
		tres=50
		x=numpy.linspace(0,R,xres)
		y=numpy.linspace(0,T,tres)
		X,Y=numpy.meshgrid(x,y)
		Z=numpy.zeros(shape=(xres,tres))
		for i in range (0,tres):
			for j in range (0,xres):
				Z[j,i]=wom_func.velocity_term(x[j],y[i])
			

		self.tck=bisplrep(X,Y,Z,s=0)
		

	def velocity(self,r,t):
		return bisplev(r,t,self.tck)

		

		
		
			
			
			
		
