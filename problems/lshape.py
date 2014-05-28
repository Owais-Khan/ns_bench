from problembase import *
from bc_ic import bc, ics
from problembase import *
from numpy import array,average

master=MPI.process_number()==0


#Fast cpp code to compute magnitude of a vector
cpp_code = """                                                                                             namespace dolfin {
	void dabla(dolfin::GenericVector& u0,dolfin::GenericVector& u1,dolfin::GenericVector& u2, dolfin::GenericVector& h,dolfin::Constant& dt, dolfin::GenericVector& cfl) {                                           	
	for (unsigned int i=0; i < cfl.size(); i++) {                                                                	cfl.setitem(i, (pow( (pow(u0[i],2) + pow(u1[i],2) +pow(u2[i],2) ) ,0.5 ))/h[i] ); 
                  }}}"""


class Problem(ProblemBase):
	def __init__(self,options):
		#Load classes
		ProblemBase.__init__(self, options)
		self.options=options
		
		# Load mesh
		self.mesh_name="data/"+self.options["mesh_name"]+".xml.gz"
		self.mesh = Mesh(self.mesh_name)

		# The body force term
		self.f = self.uConstant((0, 0, 0))

		#Viscosity
		self.nu=self.options["viscosity"]

		#Total time
		self.N_cycles=self.options["cycles"]
		self.N_timesteps=self.options["timesteps"]*self.N_cycles
		self.T=self.options["period"]*self.N_cycles
		self.dt=float(self.T)/self.N_timesteps
		self.t=0.
		self.options["current_cycle"]=0
 		self.options["file_name"]=self.options["case_name"]+"_curcyc_"+str(self.options["current_cycle"])+"_"
 
		#Save peak systolic soultion
		#by checking when Qin_new<Qin_old
		self.options["save_systole"]=True
		self.Qin_old=0 
		self.Qin_new=0
		
		#current time and cycle if restart
		if self.options["restart"]==True:
			#current time
			self.t=self.dt*float(self.options["restart_time"])+self.dt
			#current cycle
			self.options["current_cycle"]=int((((self.t/self.T)*self.N_timesteps)/self.options["timesteps"]))
			#Filename to add to store files
			self.options["file_name"]  =self.options["case_name"]+"_curcyc_"+str(self.options["current_cycle"])+"_"

		#Probe Points
		infile=open("./data/"+self.options["mesh_name"].split(".xml.gz")[0],"r")
		self.probe_points=[]
		for line in infile:
			line=line.split()
			self.probe_points.append([float(line[0]),float(line[1]),float(line[2])])
		infile.close()
		

                #Compute data to print 
                self.vol         = MPI.sum(assemble(Constant(1)*dx, mesh = self.mesh))
                self.cell_dia= [Cell(self.mesh,i).diameter() for i in range (0,self.mesh.num_cells())]
		self.avg_cell_dia=average(self.cell_dia)
		self.num_cells   = MPI.sum(self.mesh.num_cells())
		self.hmin        = MPI.min(self.mesh.hmin())
		self.hmax        = MPI.max(self.mesh.hmax())
                if master:
                        print "-"*100
                        print "Mesh Name:              ", self.mesh_name
                        print "Mesh Volume:            ", self.vol
                        print "No of cells:            ", self.num_cells
                        print "Min cell diameter:      ", self.hmin
                        print "Max cell diameter:      ", self.hmax
                        print "Average cell diameter:  ",self.avg_cell_dia
                        print "-"*100
                        print "Viscosity:              ", self.nu
                        print "Period:                 ", self.T
                        print "Time step Size:         ", self.dt
                        print "No of time steps:       ", self.N_timesteps
                        print "Current Time:           ", self.t
                        print "-"*100
                        for i in range (0,len(self.probe_points)):
                                print "Point ",i,":            ",self.probe_points[i]
                        print "-"*100







		#Calculate radius, centers, areas
		self.bc_in_ids =self.options["bc_in_ids"]
		self.bc_out_ids=self.options["bc_out_ids"]
		self.bc_in_data =self.options["bc_in_data"]
		self.bc_in_Qmean=self.options["bc_in_Qmean"]
		self.bc_func=bc.bc(self)
		self.radius, self.area, self.center=self.bc_func.geometry("inlets")
		self.radius_out, self.area_out, self.center_out=self.bc_func.geometry("outlets")
	def boundary_conditions(self, V, Q, t):		
		#Define boundary conditions
		g_noslip=self.uConstant((0,0,0))
		bc_noslip   =[DirichletBC(V,g,0) for g in g_noslip]
		bc_inflow =[]
		bc_outflow=[]
		#Inflow conditions
		for i in range(0,len(self.bc_in_ids)):
			bc_inflow.append(self.bc_func.profile(V,i,self.options["bc_in_profile"]))
		if self.options["bc_out_type"]=="zero_pressure":
			for i in range(0,len(self.bc_out_ids)):
				bc_outflow.append(self.bc_func.zero_pressure(Q,i))
		
		if self.options["bc_out_type"]=="3ewk":
			bc_outflow=self.bc_func.threeEWK(V,Q)
	
	
                #collect boundary conditions
                if len(self.bc_in_ids)==1:
                        bc_u = zip(bc_inflow[0], bc_noslip) # Important: inflow before noslip
                if len(self.bc_in_ids)==2:
                        bc_u = zip(bc_inflow[0],bc_inflow[1],bc_noslip)
                bc_p = [bc_outflow]
                return bc_u + bc_p

	

	def initial_conditions(self,V,Q):
		if self.options["restart"]:
			ics_func=ics.InitialConditions(V,Q,self.options)
			self.ics=ics_func.initial_conditions(self)
			return self.ics	
		else:
			return self.uConstant((0, 0, 0)) + [Constant(0)]

	#This function updates the current timestep
	def update(self, t, u, p):
	        self.t = t
		self.u=u
		self.p=p
		
		#for 3ewk bc
		if self.options["bc_out_type"]=="3ewk":
			for g in self.bc_func.g_outflow:
				g.update()

	

	#This function prints the output
	def functional(self,t,u,p):
		n   = tetrahedron.n
		flux_in =[]
		flux_out=[]
		pressure_in=[]
		pressure_out=[]
		try:    fd=self.mesh.domains().facet_domains()
		except: fd = MeshFunction("size_t", self.mesh, 2, self.mesh.domains())
		
		
		#Flux going in
		self.Qin_old=abs(self.Qin_new) #keep track of flux at previous ts
		for j in self.bc_in_ids:
			flux_form=Constant(0)*dx
			for i, ui in enumerate(u):
				flux_form+=inner(ui,n[i])*ds[fd](j)
			flux_in.append(assemble(flux_form))
		self.Qin_new=abs(flux_in[0]) #keep track of flux at new ts

		#Flux comming out
		for j in self.bc_out_ids:
			flux_form=Constant(0)*dx
			for i, ui in enumerate(u):
				flux_form+=inner(ui,n[i])*ds[fd](j)
			flux_out.append(assemble(flux_form))

		#Pressure at the inlets
		for j in self.bc_in_ids:
			pressure_in.append(assemble(p*ds[fd](j))/self.area[self.bc_in_ids.index(j)])
		#Pressure at the outlets
		for j in self.bc_out_ids:
			pressure_out.append(assemble(p*ds[fd](j))/self.area_out[self.bc_out_ids.index(j)])

		#Divergence of Velocity Field
                #M =  div(u)*div(u)*dx()
                #r = sqrt(assemble(M, mesh=self.mesh, form_compiler_parameters={"representation": "quadrature"}))
#		
#		#CFL Number
#		u0v=u[0].vector(); u1v=u[1].vector(); u2v=u[2].vector() 
#		func = getattr(compile_extension_module(cpp_code, "la"),  "dabla")
#		func(u0v,u1v,u2v,self.hv,Constant(self.dt),self.cflv)
#		max_CFL=self.cflv.max()
#		min_CFL=self.cflv.min()
#		avg_CFL=self.cflv.sum()/self.cflv.size()"""
		
		#Print the fluxes
		if master:

                        #Divergence
                        #print "_"*50
                        #print "Divergence of Velocity Field is ", r
#
#			#CFL 
#			print "_"*50
##			print "Maximum CFL Number in domain is ", max_CFL
#			print "Minimum CFL Number in domain is ", min_CFL
#			print "Average CFL Number in domain is ", avg_CFL"""

			
			#Flux In/Out
			print "_"*50
                        print "Flux Error is: ", 100.*(abs(sum(flux_in))-abs(sum(flux_out)))/abs(sum(flux_in)),"%"
			for i in range(0,len(self.bc_in_ids)):
				print "In   Flux     at Id ", self.bc_in_ids[i],": ",flux_in[i]
			for i in range(0,len(self.bc_out_ids)): 
				print "Out  Flux     at Id ",self.bc_out_ids[i],": ",flux_out[i]
		 	#Velocity In/out
			print "_"*50
			for i in range (0,len(self.bc_in_ids)):
				print  "In  Velocity  at Id ", self.bc_in_ids[i],": ", flux_in[i]/self.area[i]
			for i in range (0,len(self.bc_out_ids)):
				print  "Out Velocity  at Id ", self.bc_out_ids[i],": ", flux_out[i]/self.area_out[i]
			#Pressure In/out
                        print "_"*50
                        for i in range (0,len(self.bc_in_ids)):
                                print  "In  Pressure  at Id ", self.bc_in_ids[i],": ", pressure_in [i]
                        for i in range (0,len(self.bc_out_ids)):
                                print  "Out Pressure  at Id ", self.bc_out_ids[i],": ",pressure_out[i]
			print "_"*50

		#Sampe Points
		for i in range(0,len(self.probe_points)):
                	try:
                        	#try to see if it works
				if self.options["segregated"]:
                        		u0_temp=u[0](array(self.probe_points[i]))
                        		u1_temp=u[1](array(self.probe_points[i]))
                        		u2_temp=u[2](array(self.probe_points[i]))
                        		p_temp =p   (array(self.probe_points[i]))
                        		#if works, print stuff
					print "Sample Vel & Pres at Point ",i," and TS ",self.t/self.T*self.N_timesteps-1," and t=",t," is: ",u0_temp,u1_temp,u2_temp,p_temp
				else:
					u_temp=u(array(self.probe_points[i]))
					p_temp=p(array(self.probe_points[i]))
					print "Sample Vel & Pres at Point ",i," and TS ",self.t/self.T*self.N_timesteps-1," and t=",t,"is: ",u_temp[0],u_temp[1],u_temp[2],p_temp	
                	except:
                        	continue
				
			
