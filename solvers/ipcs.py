__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2008-02-01"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Anders Logg, 2008-2010.

from solverbase import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."

    def __init__(self, options):
        SolverBase.__init__(self, options)

    def solve(self, problem):

        if str(problem)=="Aneurysm":
            pc = "jacobi"
        else:
            pc = "ilu"

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = problem.timestep(problem)

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", 2)
        Q = FunctionSpace(mesh, "CG", 1)
        DG = FunctionSpace(mesh, "DG", 0)

        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)

        # Remove boundary stress term is problem is periodic
        #if is_periodic(bcp):
        #    beta = Constant(0)
        #else:
        beta = Constant(1)

        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        u0 = interpolate(u0, V)
        u1 = Function(V)
        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)
        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        # Tentative velocity step
        U = 0.5*(u0 + u)
        F1 = (1/k)*inner(v, u - u0)*dx + inner(v, grad(u0)*u0)*dx \
            + inner(epsilon(v), sigma(U, p0, nu))*dx \
            + inner(v, p0*n)*ds - beta*nu*inner(grad(U).T*n, v)*ds \
            - inner(v, f)*dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Pressure correction
        a2 = inner(grad(q), grad(p))*dx
        L2 = inner(grad(q), grad(p0))*dx - (1/k)*q*div(u1)*dx

        # Velocity correction
        a3 = inner(v, u)*dx
        L3 = inner(v, u1)*dx - k*inner(v, grad(p1 - p0))*dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Time loop
        self.start_timing()
        for t in t_range:

            # Get boundary conditions
            bcu, bcp = problem.boundary_conditions(V, Q, t)

            # Compute tentative velocity step
            b = assemble(L1)
            [bc.apply(A1, b) for bc in bcu]
            solve(A1, u1.vector(), b, "gmres", "ilu")

            # Pressure correction
            b = assemble(L2)
            #if len(bcp) == 0 or is_periodic(bcp): normalize(b)
            [bc.apply(A2, b) for bc in bcp]
            #if is_periodic(bcp):
            #    solve(A2, p1.vector(), b)
            #else:
            solve(A2, p1.vector(), b, 'gmres', 'hypre_amg')
            #if len(bcp) == 0 or is_periodic(bcp): normalize(p1.vector())

            # Velocity correction
            b = assemble(L3)
            [bc.apply(A3, b) for bc in bcu]
            solve(A3, u1.vector(), b, "gmres", pc)

            # Update
            self.update(problem, t, u1, p1)
            u0.assign(u1)
            p0.assign(p1)

        return u1, p1

    def __str__(self):
        return "IPCS"
