from dolfin import *
from scipy.special import jv
import os
import time
import numpy
import math
from numpy import array
from womersley_code import interpolate_wom
from problems import problembase

master=MPI.process_number()==0

class bc():
	def __init__(self,problem):
		self.problem=problem
		self.problembase=problembase.ProblemBase(problem.options)
		self.mesh=problem.mesh
		self.mesh_dim=self.mesh.topology().dim()
		self.options=self.problem.options
	def geometry(self,condition):
		#Inlet/Oulet of the Geometry
		#Number of Inlets and Outlets
		N_in       =len(self.problem.bc_in_ids)
		N_out      =len(self.problem.bc_out_ids)
		#Inlet and Outlet Boundary Ids
		bc_in_ids  =self.problem.bc_in_ids
		bc_out_ids =self.problem.bc_out_ids
		#Inlet and Outlet Radii and Areas
		self.radius_in  =numpy.zeros(N_in)
		self.radius_out =numpy.zeros(N_out)
		self.Area_in    =numpy.zeros(N_in)
		self.Area_out   =numpy.zeros(N_out)
		self.center_in  =numpy.zeros([N_in,3]) #center coordinates for 3 faces
		self.center_out =numpy.zeros([N_out,3])

		#Calculate the centers and radii of inlets and outlets
		VV = VectorFunctionSpace(self.mesh, "CG",1) #function space
		p = project(Expression(("x[0]", "x[1]", "x[2]")), VV)
		
		try:    fd=self.mesh.domains().facet_domains()
		except: fd = MeshFunction("size_t", self.mesh, 2, self.mesh.domains())
			
		if master: 
			print "\n"		
		#For inlets 
        	for i in range (0,N_in):
			self.Area_in[i]   = assemble(Constant(1)*ds[fd](bc_in_ids[i]),mesh=self.mesh)
			self.radius_in[i] = sqrt(self.Area_in[i]/DOLFIN_PI) 
			for j in range (0,3):
				self.center_in[i,j]=assemble(p[j]*ds[fd](bc_in_ids[i]))/self.Area_in[i]
			if master: 
				print "Inlet Radius for Id ",bc_in_ids[i]," is: ",self.radius_in[i]
				print "Inlet Area   for Id ",bc_in_ids[i]," is: ",self.Area_in[i]
		if master:
			print "\n"

		#For Outlets
	        for i in range (0,N_out):
                        self.Area_out[i]   = assemble(Constant(1)*ds[fd](bc_out_ids[i]),mesh=self.mesh)
                        self.radius_out[i] = sqrt(self.Area_out[i]/DOLFIN_PI)
                        for j in range (0,3):
		                self.center_out[i,j]=assemble(p[j]*ds[fd](bc_out_ids[i]))/self.Area_out[i]
		        if master:
				print "Outlet Radius for Id ",bc_out_ids[i]," is: ",self.radius_out[i]
        	                print "Outlet Area   for Id ",bc_out_ids[i]," is: ",self.Area_out[i]

		if condition=="inlets":
			return self.radius_in, self.Area_in, self.center_in
		if condition=="outlets":
			return self.radius_out,self.Area_out, self.center_out
	def no_slip(self,V,i):
                g_noslip=self.problembase.uConstant((0,0,0))
                bc=[DirichletBC(V,g,self.problem.bc_out_ids[i]) for g in g_noslip]
		return bc

	def zero_pressure(self,Q,face):
		bc = DirichletBC(Q,0,self.problem.bc_out_ids[face])
		return bc

	def profile(self,V,face,profile):
		self.V=V
		self.face=face
		if self.options["bc_in_profile"]=="womersley":
			self.inflow=self.problem.bc_in_data[face]
		if self.options["bc_in_profile"]=="parabolic":
			self.inflow=self.problem.bc_in_Qmean[face]
			
		class InflowComp(Expression):
			def __init__(self, V, bc, component,inflow,face):
		        	self.data = InflowData(V, bc,face,inflow)
			        self.component = component
			def eval_cell(self, values, x, ufc_cell):
				values[0] = self.data(x, ufc_cell)[self.component]
		class InflowVec(Expression):
			def __init__(self, V,bc,inflow,face):
        			self.data = InflowData(V,bc,face,inflow)
    			def eval_cell(self, values, x, ufc_cell):
        			values[:] = self.data(x, ufc_cell)
    			def value_shape(self):
        			return 3,


		class InflowData(object):
			def __init__(self,V,bc,face,inflow):
				self.face=face
				self.bc=bc
				self.center=bc.problem.center
				self.inflow=inflow
				self.mesh=V.mesh()
			       #Call the womersley class
				if profile=="womersley":
			                self.wom_prof=interpolate_wom(self.bc,self.face,self.inflow)
			def __call__(self, x, ufc_cell):
		                if profile=="womersley":
					r =sqrt( (self.center[self.face,0] - x[0])**2 + (self.center[self.face,1] - x[1])**2 + (self.center[self.face,2] - x[2])**2) 
			                val=self.wom_prof.velocity(r,self.bc.problem.t) #v(r,t)
				if profile=="plug":
					val=self.inflow/self.bc.problem.area[face]
				if profile=="parabolic":
					r =sqrt( (self.center[self.face,0] - x[0])**2 + (self.center[self.face,1] - x[1])**2 + (self.center[self.face,2] - x[2])**2) 
					val=((2.*self.inflow)*(self.bc.problem.radius[face]**2.-r**2.))/(DOLFIN_PI*self.bc.problem.radius[face]**4.)
                		cell = Cell(self.mesh, ufc_cell.index)
		                n = cell.normal(ufc_cell.local_facet)
                		return [-n.x()*val, -n.y()*val, -n.z()*val]

	        #create the inflow boundary condition
		if self.options["segregated"]==True:
			self.g_inflow=[InflowComp(V,self,d,self.inflow,face=self.face) for d in range(3)]
                else:
			self.g_inflow=[InflowVec (V,self,self.inflow,face=self.face)]
	
		bc=[DirichletBC(V,g,self.problem.bc_in_ids[self.face]) for g in self.g_inflow] 
		return bc
			


	
	def threeEWK(self,V,Q):
		class OutflowBoundaryValue(Expression):
			def init(self,bc,side,dt):
				self.problem=bc.problem
				self.bc=bc
				self.ids=self.problem.bc_out_ids
				self.current_t=-1
				self.ids=self.problem.bc_out_ids
				self.side=side
				self.R=self.bc.p2[self.side]
				self.dt=dt
				self.fd=self.problem.mesh.domains().facet_domains()

				#Assign Arrays
				self.p2   =self.bc.p2
				self.Qn   =self.bc.Qn
				self.Res  =self.bc.Res
				self.Z    =self.bc.Z
				self.C    =self.bc.C

			def eval(self, values, x):
                                values[0] = self.Qn[self.side]*self.Z[self.side] + self.p2[self.side] + (self.dt/self.C[self.side])*(self.Qn[self.side] - self.p2[self.side]/self.Res[self.side])
				
			def update(self):
				t = self.problem.t
				if t > self.current_t:
					self.R = self.resistance(t)
					self.current_t = t
					
			def resistance(self,t):
				flux = 0.0
				u = self.problem.u
				n   = tetrahedron.n
				flux_form = Constant(0)*dx
				for i, ui in enumerate(u):
					flux_form += inner(ui,n[i])*ds[self.fd](self.ids[self.side])
				self.Qn[self.side] = assemble(flux_form)
				c = self.C[self.side]
				p1 = self.Qn[self.side]*self.Z[self.side] + self.p2[self.side] + (self.dt/c)*(self.Qn[self.side] - self.p2[self.side]/self.Res[self.side])
				if master:
					print "Res coefficient at     side ", self.ids[self.side], " is ", self.Res[self.side]
				self.p2[self.side] = p1
				return p1


                bc_p=[]
                self.g_outflow = []
                ids=self.problem.bc_out_ids
                k=0;l=0

		#Assign Arrays (multiple outlets)
		self.p2   =numpy.zeros(len(ids))
                self.Qn   =numpy.zeros(len(ids))
                self.Res  =numpy.zeros(len(ids))
                self.Z    =numpy.zeros(len(ids))
		self.R_Z  =numpy.zeros(len(ids))
		self.C    =numpy.zeros(len(ids))
		
		#Parameters for MCA
		self.rho        = 1.06      #g/cm3
		self.wave_speed = 1000      #cm/s

		#compute resistences and compliances for individual outlets
		l=0
		while l<len(ids):
			R_out=self.radius_out[l]
                        self.R_Z[l]= -379097183*R_out**5 +  717408271*R_out**4 - 446219772.3*R_out**3 +123153964.2*R_out**2 - 16273155.11*R_out + 895900.6523 #g/(cm4.s)
                        self.C[l]  = 3.988103663*10**-3*R_out**4 - 5.530766856*10**-3*R_out**3 + 2.074706034*10**-3*R_out**2 - 1.462955192*10**-4*R_out + 4.606480388*10**-6 #(cm4.s2)/g

			l=l+1



                #Project old values if simulation is being restarted
                if self.options["restart"]==True:
                        #Loading Pressure
                        ics=self.problem.ics
                        if self.options["tOrder"]==1:
                                p_prev=Function(Q,ics[-1])
                        if self.options["tOrder"]>1 :
                                #Compute Previous Pressure
                                p_prev=Function(Q,ics[0][-1])
                                #Compute Previous Flow rate
                                u=as_vector([Function(V,ics[0][0]),Function(V,ics[0][1]),Function(V,ics[0][2])])
                        #Loading pressure and flowrae                                   
                        fd=self.problem.mesh.domains().facet_domains()
                        n=tetrahedron.n
                        while l<len(ids):
                                #Pressure
                                self.p2[l]=assemble(p_prev*ds[fd](ids[l]))/self.problem.area_out[l]
                                #Flow Rate
                                flux_form=Constant(0)*dx
                                for i,ui in enumerate(u):
                                        flux_form+=inner(ui,n[i])*ds[fd](ids[l])
                                self.Qn[l]=assemble(flux_form)
                                l+=1


                #Compute Resistance and loop over the outlet ids
                l=0
                while l<len(ids):
                        self.Z[l]  = (self.rho*self.wave_speed)/self.Area_out[l]
                        self.Res[l]=self.R_Z[l]-self.Z[l]
                        l+=1

                while k<len(ids):
                        g = OutflowBoundaryValue()
                        g.init(self, k, self.problem.dt)
                        bc = DirichletBC(Q, g, ids[k])
                        g.current_t = 0
                        bc_p.append(bc)
                        self.g_outflow.append(g)
                        k+=1
                return bc_p



















