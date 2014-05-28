__author__ = "Harish Narayanan <harish@simula.no>"
__date__ = "2009-01-20"
__copyright__ = "Copyright (C) 2009-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kristian Valen-Sendstad, 2009.
# Modified by Anders Logg, 2010.

from problembase import *
from numpy import array
from math import pi, e

# Problem definition
class Problem(ProblemBase):
    "3D test problem with known analytical solution."

    def __init__(self, options):
        ProblemBase.__init__(self, options)

        # We start with a UnitCube and modify it to get the mesh we
        # want: (-1, 1) x (-1, 1) x (-1, 1)

        mesh_sizes = [5, 8, 11, 16, 23, 32]
        level = options["refinement_level"]
        N = int(mesh_sizes[level])

        self.mesh = UnitCube(N, N, N)
        self.scale  = 2*(self.mesh.coordinates() - 0.5)
        self.mesh.coordinates()[:, :] = self.scale

        # The body force term
        self.f = self.uConstant((0, 0, 0))

        # Set the kinematic viscosity (nu = eta/rho)
        self.nu = 1.0

        # Characteristic velocity in the domain (used to determine timestep)
        self.U = 1.0

        # Set final time
        self.T = 0.5

	#Time
	self.dt, self.t, self.t_range =self.timestep(self)

        # The analytical solution
        # Velocity
        self.analytical_u = \
            ('-((a*(pow(E,a*x[2])*cos(a*x[0] + d*x[1]) + pow(E,a*x[0])*sin(a*x[1] + d*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
             '-((a*(pow(E,a*x[0])*cos(a*x[1] + d*x[2]) + pow(E,a*x[1])*sin(d*x[0] + a*x[2])))/pow(E,pow(d,2)*t*etabyrho))',
             '-((a*(pow(E,a*x[1])*cos(d*x[0] + a*x[2]) + pow(E,a*x[2])*sin(a*x[0] + d*x[1])))/pow(E,pow(d,2)*t*etabyrho))')
        # Pressure
        self.analytical_p = \
            ('-(rho/2.0)*(pow(a,2)*(pow(E,2*a*x[0]) + pow(E,2*a*x[1]) + pow(E,2*a*x[2]) + 2*pow(E,a*(x[1] + x[2]))*cos(d*x[0] + a*x[2])*sin(a*x[0] + d*x[1]) + 2*pow(E,a*(x[0] + x[1]))*cos(a*x[1] + d*x[2])*sin(d*x[0] + a*x[2]) + 2*pow(E,a*(x[0] + x[2]))*cos(a*x[0] + d*x[1])*sin(a*x[1] + d*x[2])))/(pow(E,pow(d,2)*t*etabyrho))')

        # Common parameters pertinent to the functional forms above
        self.u_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e,             'etabyrho': 1.0, 't': 0.0}
        self.p_params = {'a': pi/4.0, 'd': pi/2.0, 'E': e, 'rho': 1.0, 'etabyrho': 1.0, 't': 0.0}

    def initial_conditions(self, V, Q):

        # Use analytical solutions at t = 0 as initial values
        self.exact_u = self.uExpr(self.analytical_u, degree=3, **self.u_params)
        self.exact_p = self.uExpr(self.analytical_p, degree=3, **self.p_params)

        return self.exact_u + self.exact_p

    def boundary_conditions(self, V, Q, t):
        for expr in self.exact_u + self.exact_p:
            expr.t = t

        bc0 = [DirichletBC(V, expr, DomainBoundary()) for expr in self.exact_u]

        bcu   = zip(bc0)
        bcp   = [()]

        return bcu + bcp

    def update(self, t, u, p):
        for expr in self.exact_u + self.exact_p:
            expr.t = t

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

    def reference(self, t):
        return 0.0

    def __str__(self):
        return "Beltrami"
