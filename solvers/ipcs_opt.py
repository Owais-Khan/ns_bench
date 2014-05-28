from __future__ import division

__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2011-11-11"
__copyright__ = "Copyright (C) 2011 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from solverbase import *
from rhsgenerator import *

class Solver(SolverBase):
    "Incremental pressure-correction scheme."
    set_log_active(False)
    def __init__(self, options):
        SolverBase.__init__(self, options)
        self.segregated = options['segregated']

    def solve(self, problem, restart=None):

        # Get problem parameters
        mesh = problem.mesh
        dt, t, t_range = self.select_timestep(problem)
        if restart:
            dt, t, t_range = restart.select_timestep(dt, problem.T)

        # Define function spaces
        if self.segregated:
            V = FunctionSpace(mesh, "CG", 1)
        else:
            V = VectorFunctionSpace(mesh, "CG", 1)
        Q = FunctionSpace(mesh, "CG", 1)

        # Get initial and boundary conditions
        ics = problem.initial_conditions(V, Q)
        if restart:
            u0 = restart.u(t, V)
            p0 = restart.p(t, Q)
        else:
            u0 = [interpolate(_, V) for _ in ics[:-1]]
            p0 = interpolate(ics[-1], Q)
        bcs = problem.boundary_conditions(V, Q, t)
        bcu, bcp = bcs[:-1], bcs[-1]

        # Remove boundary stress term if problem is periodic
        #beta = 0 if is_periodic(bcp) else 1
	beta=1
        # Test and trial functions
        v = TestFunction(V)
        q = TestFunction(Q)
        u = TrialFunction(V)
        p = TrialFunction(Q)

        # Functions
        dims = range(len(u0))
        u1 = [Function(V) for d in dims]
        p1 = Function(Q)

        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f
        n  = FacetNormal(mesh)

        if not self.segregated:
            # To avoid indexing in non-segregated forms
            u0_ = u0[0]
            u1_ = u1[0]
            f_  = f[0]

        # Tentative velocity step
        M  = assemble(inner(v, u) * dx)
        K1 = assemble((1/k) * inner(v, u) * dx)
        if self.segregated:
            K2 = assemble(0.5 * inner(grad(v), nu*grad(u)) * dx)
            a_conv = -v * sum(u0[r]*u.dx(r) for r in dims) * dx
            Kconv = Matrix() # assembled from a_conv in the time loop from

            A_u_tent = []
            rhs_u_tent = []
            for d in dims:
                A_u_tent.append(K1+K2) # Separate matrices, because they may have different BCs
                K3 = assemble(-v*p*n[d]*ds + v.dx(d)*p*dx)

                rhs = RhsGenerator(V)
                rhs += K1, u0[d]
                rhs -= K2, u0[d]
                rhs += K3, p0
                rhs += M, f[d]
                rhs += Kconv, u0[d]

                rhs_u_tent.append(rhs)
        else:
            K2 = assemble(inner(epsilon(v), nu*epsilon(u)) * dx
                          - 0.5 * beta * nu * inner(v, grad(u).T*n) * ds)
            A = K1+K2
            K3 = assemble(-inner(v, p*n)*ds + div(v)*p*dx)

            rhs = RhsGenerator(V)
            rhs += K1, u0_
            rhs -= K2, u0_
            rhs += K3, p0
            rhs += M, f_
            rhs += -inner(v, grad(u0_)*u0_) * dx

            A_u_tent, rhs_u_tent = [A], [rhs]

        # Pressure correction
        A_p_corr = assemble(inner(grad(q), grad(p))*dx)
        if self.segregated:
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            for r in dims:
                Ku = assemble(-(1/k)*q*u.dx(r)*dx)
                rhs_p_corr += Ku, u1[r]
        else:
            Ku = assemble(-(1/k)*q*div(u)*dx)
            rhs_p_corr = RhsGenerator(Q)
            rhs_p_corr += A_p_corr, p0
            rhs_p_corr += Ku, u1_

        # Velocity correction
        A_u_corr = [M.copy() for r in dims]
        if self.segregated:
            rhs_u_corr = []
            for r in dims:
                Kp = assemble(-k*inner(v, grad(p)[r])*dx)
                rhs = RhsGenerator(V)
                rhs += M, u1[r]
                rhs += Kp, p1
                rhs -= Kp, p0
                rhs_u_corr.append(rhs)
        else:
            Kp = assemble(-k*inner(v, grad(p))*dx)
            rhs_u_corr = RhsGenerator(V)
            rhs_u_corr += M, u1_
            rhs_u_corr += Kp, p1
            rhs_u_corr -= Kp, p0
            rhs_u_corr = [rhs_u_corr]

        # Apply BCs to matrices
        for A, bcs in zip(A_u_tent, bcu) + zip(A_u_corr, bcu):
            for bc in bcs:
                bc.apply(A)
        for bc in bcp:
            bc.apply(A_p_corr)

        # Create solvers
        if len(bcp)==0:# or is_periodic(bcp):
            solver_p_params = self.options['solver.p'] or self.options['solver.p_neumann']
        else:
            solver_p_params = self.options['solver.p'] or self.options['solver.p_dirichlet']
        solver_u_tent = [LinearSolver(*self.options['solver.u_tent']) for d in dims]
        solver_p_corr = LinearSolver(*solver_p_params)
        solver_u_corr = [LinearSolver(*self.options['solver.u_corr']) for d in dims]

        for A,S in zip(A_u_tent, solver_u_tent) \
                + [(A_p_corr, solver_p_corr)] \
                + zip(A_u_corr, solver_u_corr):
            S.set_operator(A)
            if 'preconditioner' in S.parameters:
                S.parameters['preconditioner']['reuse'] = True

        # Time loop
        self.start_timing()
        for t in t_range:
            # Get boundary conditions
            bcs = problem.boundary_conditions(V, Q, t)
            bcu, bcp = bcs[:-1], bcs[-1]
            self.timer("update & fetch bc")

            # Assemble the u0-dependent convection matrix. It is important that
            # it is assembled into the same tensor, because it is stored in rhs.
            if self.segregated:
                assemble(a_conv, tensor=Kconv, reset_sparsity=(Kconv.size(0)==0))

            # Compute tentative velocity step
            for d, S, rhs, u1_comp, bcu_comp in zip(dims, solver_u_tent, rhs_u_tent, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                self.timer("u0 construct rhs")
                iter = S.solve(u1_comp.vector(), b)
                self.timer("u0 solve (%s, %d, %d)"%(', '.join(self.options['solver.u_tent']), A.size(0), iter))

            # Pressure correction
            b = rhs_p_corr()
            if len(bcp) == 0: #or is_periodic(bcp): 
		normalize(b)
            for bc in bcp: bc.apply(b)
            self.timer("p1 construct rhs")
            iter = solver_p_corr.solve(p1.vector(), b)
            if len(bcp) == 0: #or is_periodic(bcp): 
		normalize(p1.vector())
            self.timer("p1 solve (%s, %d, %d)"%(', '.join(solver_p_params), A_p_corr.size(0), iter))

            # Velocity correction
            for S, rhs, u1_comp, bcu_comp in zip(solver_u_corr, rhs_u_corr, u1, bcu):
                b = rhs()
                for bc in bcu_comp: bc.apply(b)
                self.timer("u1 construct rhs")
                iter = S.solve(u1_comp.vector(), b)
                self.timer("u1 solve (%s, %d, %d)"%(', '.join(self.options['solver.u_corr']), A.size(0),iter))

            # Update
            self.update(problem, t, u1, p1)
            for r in dims: u0[r].assign(u1[r])
            p0.assign(p1)

        return as_object(u1), p1

    def __str__(self):
        name = "IPCS_opt"
        if self.segregated:
            name += "_seg"
        return name
