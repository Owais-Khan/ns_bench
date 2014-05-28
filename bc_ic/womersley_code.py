import math
import numpy
from scipy.special import jv
from scipy import *
from numpy import array
import os
from time import time
from scipy.interpolate import splrep, splev, interp2d, interpolate, bisplrep, bisplev
from dolfin import *
master = MPI.process_number() == 0

class flow_waveform():
	def __init__(self,bc,face,FC):
		#Input parameters
		self.problem=bc.problem
		self.N_cycles=self.problem.N_cycles
		self.N_timesteps=self.problem.N_timesteps
		self.nu=self.problem.nu
		self.R=self.problem.radius[face]
		self.T0=self.problem.T
		self.dt=self.problem.dt
		self.options=self.problem.options
		
		#Load the Fourier Coefficients
		infile_FC=open("./data/"+FC)
		self.an=[];self.bn=[]
		for line in infile_FC:
			self.an.append(float(line.split()[0]))
			self.bn.append(float(line.split()[1]))
			
		self.omega=((2.*DOLFIN_PI)/self.T0)*self.N_cycles #omega
		self.Q_mean=self.problem.bc_in_Qmean[face] #Mean flow rate

		self.womn=self.R*(self.omega/self.nu)**0.5 #womersley number
		if master:
			print "Womersley Number is: :",self.womn

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

		
		#Directory to store the velocity profile
		self.wom_dir="WomProf_"+self.options["case_name"]+"_womn"+str(self.womn)#Qmean"+str(self.an[0]*self.Q_mean)+"_R"+str(self.R)+"_nu"+str(self.nu)+"_ts"+str(self.N_timesteps/self.N_cycles)+"_cycles"+str(self.N_cycles)		
		self.n=0

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
		self.problem=bc.problem
		wom_prof=flow_waveform(bc,face,FC)
		xres=50
		timesteps=int(self.problem.N_timesteps/self.problem.N_cycles)
		self.t=self.problem.t
		
		wom_dir=wom_prof.wom_dir

		#define the spatial and temporal timesteps
		self.fr=numpy.linspace(0,self.problem.radius[face],xres,endpoint=True)
		
		#open a file to store the value
		if wom_dir in os.listdir("./data/"):abcd=1
		else:
			self.infile=open("./data/"+wom_dir,'w')
			for i in range (0,timesteps):
				for j in range (0,len(self.fr)):
					self.infile.write(str(wom_prof.velocity_term(self.fr[j],self.t))+" ")
				if i<timesteps-1: self.infile.write("\n")
				self.t=self.t+self.problem.dt
		        self.infile.close()

		#reset values
		self.infile=open("./data/"+wom_dir,'r')
		self.t=-1
		
		#if t starts somewhere in the middle
		if self.problem.options["restart"]: 
			self.time_step=self.problem.options["restart_time"]%self.problem.options["timesteps"]	
			for i in range (0,self.time_step+1): self.infile.readline()

	def velocity(self,r,t):
		if self.t!=t:
			if self.problem.options["current_timestep"]%self.problem.options["timesteps"]==self.problem.options["timesteps"]-1:
				self.infile.seek(0,0)	
			line=self.infile.readline()
			line=line.split()
			val=[]
			for i in line: val.append(float(i))
			self.tck=splrep(self.fr,val)
			self.t=t
			return splev(r,self.tck)
		else:
			return splev(r,self.tck)

		

		
		
			
			
			
		
