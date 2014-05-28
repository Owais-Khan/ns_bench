__author__ = "Muhammad Owais Khan <owais.khan@alumni.utoronto.ca>"
__date__ = "2013-12-11"
__license__  = "GNU GPL version 3 or any later version"


from problembase import *
from numpy import array
master=MPI.process_number()==0


#Define the domain boundaries
class PeriodicDomain(SubDomain):
    
    def inside(self, x, on_boundary):
        # return True if on left or bottom boundary AND NOT on one of the two corners (0, 2) and (2, 0)
        return bool((near(x[0], 0) or near(x[1], 0)) and
              (not ((near(x[0], 0) and near(x[1], 2)) or
                    (near(x[0], 2) and near(x[1], 0)))) and on_boundary)

    def map(self, x, y):
        if near(x[0], 2) and near(x[1], 2):
            y[0] = x[0] - 2.0
            y[1] = x[1] - 2.0
        elif near(x[0], 2):
            y[0] = x[0] - 2.0
            y[1] = x[1]
        else:
            y[0] = x[0]
            y[1] = x[1] - 2.0


# Problem definition
class Problem(ProblemBase):
    "2D Taylor Green"

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = RectangleMesh(0,0,2,2,N,N,"crossed")

        # Create right-hand side function with pressure gradient as body force
        self.f = self.uConstant((0, 0))

        # Parameters 
        self.nu       = 0.01   #viscosity
	self.T        = 1.     #period
	self.timesteps= 1000    #timesteps
	self.dt       = self.T/self.timesteps
	self.t        = 0.0	

	#define analytical solution
	self.analytical_u=('-sin(pi*x[1])*cos(pi*x[0])*exp(-2.*pi*pi*nu*t)',
			   ' sin(pi*x[0])*cos(pi*x[1])*exp(-2.*pi*pi*nu*t)')
	self.analytical_p=('-(cos(2*pi*x[0])+cos(2*pi*x[1]))*exp(-4.*pi*pi*nu*t)/4.')
	self.params      ={'nu': self.nu,'t':0.0}
	
	#Periodic domain
	options["constrained_domain"]=PeriodicDomain()
	if master: "Pring using Periodic Domain"
	
    def initial_conditions(self, V, Q):
        # Use analytical solutions at t = 0 as initial values
        self.exact_u = self.uExpr(self.analytical_u, degree=2, **self.params)
        self.exact_p = self.uExpr(self.analytical_p, degree=2, **self.params)
        return self.exact_u + self.exact_p

    def boundary_conditions(self, V, Q, t):

        # Create no-slip boundary condition for velocity
	bcu=[(),()]
	bcp=[()]
        return bcu + bcp

    def functional(self, t, u, p):
        if t < self.T:
            return 0.0
        else:
            for expr in self.exact_u + self.exact_p:
                expr.t = t
            if not self.options['segregated']:
                u = [u]
            error = 0
            for exact_u, calc_u in zip(self.exact_u, u):
                error += sqr(errornorm(exact_u, calc_u) / norm(exact_u, mesh=self.mesh))
            return sqrt(error/len(u))

    	return sqrt(error/len(u))



    def reference(self, t):
        return None

    def tolerance(self, problem):
        return 1e-11

    def __str__(self):
        return "Taylor Green 2d"
