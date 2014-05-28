__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from problembase import *

# Boundary value
def boundaryvalue(x):
    if x[0] > DOLFIN_EPS and x[0] < 1.0 - DOLFIN_EPS and x[1] > 1.0 - DOLFIN_EPS:
        return [1.0, 0.0]
    else:
        return [0.0, 0.0]

class BoundaryValueVec(Expression):
    def value_shape(self):
        return (2,)
    def eval(self, values, x):
        values[:] = boundaryvalue(x)

class BoundaryValueComp(Expression):
    def __init__(self, component, **kwargs):
        Expression.__init__(self, **kwargs)
        self.component = component
    def eval(self, values, x):
        values[0] = boundaryvalue(x)[self.component]

# Problem definition
class Problem(ProblemBase):
    "2D lid-driven cavity test problem with known reference value."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # Create mesh
        N = options["N"]
        self.mesh = UnitSquare(N, N)

        # Create right-hand side function
        self.f = self.uConstant((0, 0))

        # Set viscosity (Re = 1000)
        self.nu = 1.0 / 1000.0
        self.U = 1.0

        # Set end-time
        self.T = 2.5

	#dt
	self.dt=0.025

	#timestels
	self.time_steps=self.T/self.dt

	#current t
	self.t=0.0

    def initial_conditions(self, V, Q):

        u0 = self.uConstant((0, 0))
        p0 = self.uConstant(0)

        return u0 + p0

    def boundary_conditions(self, V, Q, t):

        if self.options['segregated']:
            element = FiniteElement("CG", triangle, 1)
            self.g = [BoundaryValueComp(d, element=element) for d in range(2)]
        else:
            element = VectorElement("CG", triangle, 1)
            self.g = [BoundaryValueVec(element=element)]
        bc = [DirichletBC(V, g, DomainBoundary()) for g in self.g]

        return zip(bc) + [()]

    def functional(self, t, u, p):

        # Only check final time
        if t < self.T:
            return 0
        else:
            # Compute stream function and report minimum
            psi = StreamFunction(u)
            vals  = psi.vector().array()
            vmin = MPI.min(vals.min())

            if master:
                print "Stream function has minimal value" , vmin

            return vmin

    def reference(self, t):

        # Only check final time
        if t < self.T:
            return 0.0

        return -0.061076605

    def __str__(self):
        return "Driven cavity"

def StreamFunction(u):
    "Stream function for a given 2D velocity field."

    # Fetch a scalar function (sub-)space
    try:
        V = u.function_space()
        V = V.sub(0).collapse()
    except AttributeError:
        V = u[0].function_space()

    # Check dimension
    mesh = V.mesh()
    if not mesh.topology().dim() == 2:
        error("Stream-function can only be computed in 2D.")

    # Define variational problem
    q   = TestFunction(V)
    psi = TrialFunction(V)
    a   = dot(grad(q), grad(psi))*dx
    L   = dot(q, (u[1].dx(0) - u[0].dx(1)))*dx

    # Define boundary condition
    g  = Constant(0)
    bc = DirichletBC(V, g, DomainBoundary())

    # Compute solution
    psi = Function(V)
    solve(a == L, psi, bc)

    return psi
