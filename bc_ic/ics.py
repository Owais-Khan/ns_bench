from dolfin import *
master=MPI.process_number()==0
from problems import problembase
import os
class InitialConditions():

	def __init__(self,V,Q,options):
		self.problembase=problembase.ProblemBase(options)
		self.V=V
		self.Q=Q
		self.mesh=V.mesh()
		self.options=options

	def initial_conditions(self,problem):        

		#File paths     	    
                restart_path="./results/"+self.options["case_name"]+"/check_point"
                files=sorted(os.listdir(restart_path),key=lambda p: os.path.getctime(os.path.join(restart_path,p)))

		#Create the velocity file names
		ufiles_1=[]; ufiles_2=[]
		if   self.options["segregated"]:n=3
		else:n=1
		restart_time_1=str("%.5d"%(self.options["restart_time"]))
		restart_time_2=str("%.5d"%(self.options["restart_time"]-1))
		for i in range(n):
			for j in files:
							
				if j.find("u"+str(i)+"_ts="+restart_time_1)>=0: 
					ufiles_1.append(restart_path+"/"+j)
				if j.find("u"+str(i)+"_ts="+restart_time_2)>=0:
					ufiles_2.append(restart_path+"/"+j)
		
		#create the pressure file names
		pfile_1=[]; pfile_2=[]
                for i in range(1):
                        for j in files:
                                if j.find("p_ts="+restart_time_1)>=0: pfile_1.append(restart_path+"/"+j)
                                if j.find("p_ts="+restart_time_2)>=0: pfile_2.append(restart_path+"/"+j)



                #Define a new function space            
                self.V2 = FunctionSpace(self.mesh, "CG", self.options["uOrder"]) # The new function space
                self.Q2 = FunctionSpace(self.mesh, "CG", 1)
                u = Function(self.V)
                p = Function(self.Q)


		#Apply the files to function space
		ics=[]; ics_temp=[]
		for i in range(n):
			if master: print 'Projecting ', ufiles_1[i]
			ics_temp.append(Function(self.V,ufiles_1[i]))
		if master: print 'Projecting ', pfile_1[0]
		ics_temp.append(Function(self.Q,pfile_1[0]))
		ics.append(ics_temp);ics_temp=[]
		
		if self.options["segregated"]:
	                for i in range(n):
        	                if master: print 'Projecting ', ufiles_2[i]
                	        ics_temp.append(Function(self.V,ufiles_2[i]))
	                if master: print 'Projecting ', pfile_2[0]
        	        ics_temp.append(Function(self.Q,pfile_2[0]))
                	ics.append(ics_temp)


		if self.options["tOrder"]==1: return ics[0]			
		if self.options["tOrder"]==2: return ics



					
	
