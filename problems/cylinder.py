__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2009-10-01"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2010.

from problembase import *
from numpy import array
import matplotlib.pyplot as plt
import os

# Constants related to the geometry
bmarg = 1.e-3 + DOLFIN_EPS
xmin = 0.0
xmax = 2.2
ymin = 0.0
ymax = 0.41
xcenter = 0.2
ycenter = 0.2
radius = 0.05

#Cylinder Domain
class Cylinder(SubDomain):
    def inside(self,x,on_boundary):
	dx=x[0]-xcenter
	dy=x[1]-ycenter
	r =sqrt(dx*dx+dy*dy)
	return on_boundary and r<radius+bmarg	

# Inflow boundary
class InflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] < xmin + bmarg

# No-slip boundary
class NoslipBoundary(SubDomain):
    def inside(self, x, on_boundary):
        dx = x[0] - xcenter
        dy = x[1] - ycenter
        r = sqrt(dx*dx + dy*dy)
        return on_boundary and \
               (x[1] < ymin + bmarg or x[1] > ymax - bmarg or \
                r < radius + bmarg)

# Outflow boundary
class OutflowBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and x[0] > xmax - bmarg

# Problem definition
class Problem(ProblemBase):

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Load mesh
        self.refinement_level = options["refinement_level"]
        if self.refinement_level > 5:
            raise RuntimeError, "No mesh available for refinement level %d" % self.refinement_level

        self.mesh = Mesh("data/cylinder_%d.xml.gz" % self.refinement_level)


        # Create right-hand side function
        self.f = self.uConstant((0,0))

	#Set problem parameters (Schafer and Turek)
	self.Re    = 100                        #Reynolds number (20-steady or 100-unsteady)
        self.nu    = 1.0 / 1000.0               #Kinematic viscosity
	self.rho   = 1.                         #Density
	self.mu    = self.nu*self.rho           #Dynamic viscoisty
	self.u_bar = (self.Re*self.nu)\
		    /(self.rho*2*radius)        #Mean velocity
	self.Umax  = 3./2.*self.u_bar           #Maximum velocity
	self.U     = self.Umax                  #Time stepping velocity
	self.T     = 8.0                        #End Time
#	self.dt    = 1.*self.mesh.hmin()\
#		    /self.Umax	                #dt based on CFL=0.025
#	self.Nts   = int(self.T/self.dt)        #Number of Time steps
	self.Nts   = 1000
	self.dt    = self.T/self.Nts

	#Update case name from default
	options["case_name"]=str(options["problem_name"])+"_"+str(options["solver_name"])+"_refinement_level"+\
			    str(self.refinement_level)+"_ts"+str(self.Nts)+"_cycles"+str(1)\
			    +"_uOrder"+str(options["uOrder"])
	self.case_name     =options["case_name"]
	self.casedir = os.path.join("results",self.case_name)
        if not os.path.exists(self.casedir): os.mkdir(self.casedir)

	##########Post-Processing########## 
	
	#Normal Facet and markers to compute
	#drag and lift over cylinder
	self.n       = FacetNormal(self.mesh)
	self.markers = MeshFunction("size_t", self.mesh, self.mesh.topology().dim() - 1) 
	Cylinder().mark(self.markers, 100)

	#Compute location of end of recirculation zone
	#Search between end of cylinder (0.25m) to 0.5m
	#Only valid for Re=20 not Karmen vortex street(Re=100)
	self.La_array=linspace(0.25,0.5,500)

	#Create a Tecplot Readable ASCII file:
	#Store Time, Pressure Drop, Length of Recirculation zone, drag and lift
	file_name=self.casedir+"/Cylinder_refinement_level_"+str(self.refinement_level)+"_Re"+str(self.Re)
	self.outfile=open(file_name+".dat",'w')
        self.outfile.write('TITLE = "'+file_name+'"\n')
        self.outfile.write('VARIABLES = "t", "dp", "La", "Cd", "Cl"\n')
        self.outfile.write('Zone T= "'+file_name+'", I='+str(self.Nts)+', F=POINT\n')

	
    def initial_conditions(self, V, Q):
        u0 = self.uConstant((0,0))
        p0 = self.uConstant(0)

        return self.uConstant((0, 0)) + [Constant(0)]

    def boundary_conditions(self, V, Q, t):
        # Create inflow boundary condition
        self.b0 = InflowBoundary()
        self.g0 = self.uExpr(('4*Um*(x[1]*(ymax-x[1]))/(ymax*ymax)', '0.0'),
                             Um=self.Umax, ymax=ymax, t=t)

        bc0 = [DirichletBC(V, g0, self.b0) for g0 in self.g0]

        # Create no-slip boundary condition
        self.b1 = NoslipBoundary()
        self.g1 = self.uConstant((0, 0))
        bc1     = [DirichletBC(V, g1, self.b1) for g1 in self.g1]

        # Create outflow boundary condition for pressure
        self.b2 = OutflowBoundary()
        self.g2 = Constant(0)
        bc2     = [DirichletBC(Q, self.g2, self.b2)]

        # Collect boundary conditions
        bcu = zip(bc0, bc1)
        bcp = zip(bc2)

        return bcu + bcp

    def update(self, t, u, p):
        for g0 in self.g0:
            g0.t = t

    def functional(self, t, u, p):
	#Compute Pressure difference	
        x1 = array((0.15,0.2)); x2 = array((0.25,0.2))
	dp= p(x1)-p(x2)

	#Compute the length of recirculation zone 
	#When x-velocity becomes positive in wake
	for xr in self.La_array:
		if u[0](array((xr,ymax/2.)))>=0:
			La_u=u[0](array((xr,ymax/2.)))
			break

	#Compute lift and drag
	drag,lift=self.drag_lift(t,u,p)

	#Print data
	print "The Pressure difference is:          ",p(x1)-p(x2)
	print "The Length of recirculation zone is: ",xr-0.25
	print "The Drag coefficient is            : ",drag
	print "The Lift coefficient is            : ",lift

	#Write to outfile
	self.outfile.write(str(t)+" "+str(p(x1)-p(x2))+" "+str(xr-0.25)+" "+str(drag)+" "+str(lift)+"\n")

	#Compute stream Function
	self.stream_function(t,u,p)

        return 0


    def drag_lift(self,t,u,p):
	    R = VectorFunctionSpace(self.mesh, 'R', 0)
	    c = TestFunction(R) 
	    n = FacetNormal(self.mesh)	    
	    tau = self.nu * (nabla_grad(u) + nabla_grad(u).T) - p*Identity(p.cell().d)
	    F   = assemble(dot(dot(tau, n), c)*ds(100),exterior_facet_domains=self.markers)
	    F   = (F)/(self.u_bar**2*radius)
	    return F[0],F[1]

    def stream_function(self,t,u,p):
	    # StreamFunction and vorticity
	    F = FunctionSpace(self.mesh,'CG',2)
	    w = project(curl(u),F)

	    #Define trial and test function
	    psi  = TrialFunction(F)
	    psi_v = TestFunction(F)

	    #Facet normals and variational form
	    n = FacetNormal(self.mesh)
	    a = -inner(grad(psi ),grad(psi_v))*dx
	    l = -w*psi_v*dx + (n[0]*u[1] - n[1]*u[0])*psi_v*ds

	    #Assemble and solve the equations
	    A,b = assemble_system(a,l)
	    psi  = Function(F)
	    solve(A,psi.vector(),b)

	    #Normalize the vector and reorient the array
	    normalize(psi.vector())
	    psi_min= abs(psi.vector().array().min())
	    psi.vector()[:]+=psi_min

	    #Export the solution
	    output = File(self.casedir+"/StreamFunction_refinement_level_"+str(self.refinement_level)+"_Re"+str(self.Re)+".pvd")   
	    output <<psi


    def reference(self, t):
	return None

    def tolerance(self, problem):
        return 1e-7

    def __str__(self):
        return "Cylinder"
