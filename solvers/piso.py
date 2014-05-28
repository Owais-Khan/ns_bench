__author__ = "Kristian Valen-Sendstad <kvs@simula.no>"
__date__ = "2012"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"


from solverbase import *
class Solver(SolverBase):
    "PISO, Issa 1985, implemented according to algorithm in Versteeg and Malalasekera"

    def __init__(self, options):
        assert not options['segregated']
        SolverBase.__init__(self, options)

    def solve(self, problem):

        parameters["linear_algebra_backend"] = "PETSc"
        parameters["form_compiler"]["optimize"]     = True
        parameters["form_compiler"]["cpp_optimize"] = True


        solver_p_periodic  = "gmres", "hypre_euclid"

        solver_p_dirichlet = "gmres", "ml_amg"#, "parameters['krylov_solver']['preconditioner']['reuse'] = True"
        solver_p_corr = LinearSolver(*solver_p_dirichlet)
        solver_p_corr.parameters['preconditioner']['reuse'] = True
        solver_p_corr.parameters['nonzero_initial_guess'] = True

        solver_u_prm      = "bicgstab", "hypre_euclid"
        solver_u_corr = LinearSolver(*solver_u_prm)
        solver_u_corr.parameters['preconditioner']['reuse'] = True
        solver_u_corr.parameters['nonzero_initial_guess'] = True
        # Get problem parameters
        mesh = problem.mesh
        #dt, t, t_range = problem.timestep(problem)
	dt=problem.dt
	t=problem.t
	T=problem.T

        # Define function spaces
        V = VectorFunctionSpace(mesh, "CG", self.options["uOrder"])
        Q = FunctionSpace(mesh, "CG", 1)


        # Get initial and boundary conditions
        u0, p0 = problem.initial_conditions(V, Q)
        bcu, bcp = problem.boundary_conditions(V, Q, t)


        # Remove boundary stress term if problem is periodic
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
	u0  = interpolate(u0,V)
	u2  = interpolate(u0,V)
	un  = interpolate(u0,V)
	unc = interpolate(u0,V)
	#u0 = Function(V,u0)
	#u2 = Function(V)
        #un = Function(V)
        #unc = Function(V)

        p0 = interpolate(p0, Q)
        p1 = interpolate(p0, Q)

        nu = Constant(problem.nu)
        k  = Constant(dt)
        f  = problem.f[0]
        n  = FacetNormal(mesh)


        num_dofs = u0.vector().size()
        print 'num u dofs is ', num_dofs
        num_dofs = u0.vector().size() + p0.vector().size()
        print 'tot num dofs is ', num_dofs

        if False: # Cebral implementation of SUPG parameter
            dl = ((12/sqrt(2))*(tetrahedron.volume))**.33333333
            supg =   (1/k)*(dl**2/(2*sqrt(inner(u0,u0))*dl+4*nu))
            supgII = (1/k)*(dl**2/(2*sqrt(inner(u2,u2))*dl+4*nu))
            print 'Using supg with dl**2/(2*sqrt(inner(u,u))*dl+4*nu)/k'
        else:
            supg = Constant(1)
            supgII = Constant(1)
            print 'Using supg with Constant(1)'

        # Predictor
        F_u_tent = (1/k)*inner(v, un - u0)*dx  + inner(v,grad(p0))*dx + inner(v, grad(un)*un)*dx + nu*inner(grad(v), grad(un))*dx \
            - inner(v, f)*dx + supg*k*inner(grad(v)*un, grad(un)*(un))*dx

        # C1, Pressure correction I
        a_p_corr = inner(grad(q), grad(p))*dx
        L_p_corr = inner(grad(q), grad(p0))*dx  - (1/k)*q*div(un)*dx

        # C1, Velocity correction I
        a_u_corr = inner(v, u)*dx
        L_u_corr = inner(v, un)*dx - k*inner(v, grad(p1-p0))*dx

        # C2, Pressure correction II
        a_p_corrII = inner(grad(q), grad(p))*dx
        L_p_corrII = inner(grad(q), grad(p1))*dx - (1/k)*q*div(u2)*dx

        # C2, Velocity correction I
        F_u_corrm = (1/k)*inner(v, unc - u0)*dx + inner(v, grad(p1))*dx + inner(v, grad(unc)*unc)*dx  + nu*inner(grad(v), grad(unc))*dx \
            - inner(v, f)*dx + supgII*k*inner(grad(v)*unc, grad(unc)*(unc))*dx


        # Assemble matrices
        A_p_corr = assemble(a_p_corr)
        A_u_corr = assemble(a_u_corr)
        A_p_corrII = assemble(a_p_corrII)
        print 'assemble  ok'


        J = derivative(F_u_tent, un, u)
        uproblem = NonlinearVariationalProblem(F_u_tent, un, bcs = list(bcu), J = J)
        usolver = NonlinearVariationalSolver(uproblem)
        prm = usolver.parameters
        prm["newton_solver"]["absolute_tolerance"] = 1E-6
        prm["newton_solver"]["relative_tolerance"] = 1E-6
        prm["newton_solver"]["maximum_iterations"] = 1
        prm["newton_solver"]["relaxation_parameter"] = .99
        prm["newton_solver"]["error_on_nonconvergence"] = False
        prm['linear_solver'] = 'gmres'
        prm['preconditioner'] = 'hypre_euclid'
#        prm['linear_solver'] = 'bicgstab'
#        prm['preconditioner'] = 'jacobi'
        prm['krylov_solver']['preconditioner']['reuse'] = True
        prm['krylov_solver']['absolute_tolerance'] = 1E-7
        prm['krylov_solver']['relative_tolerance'] = 1E-6
        prm['krylov_solver']['maximum_iterations'] = 20000
        prm['krylov_solver']['monitor_convergence'] = True
        prm['krylov_solver']['nonzero_initial_guess'] = True
        prm['krylov_solver']['gmres']['restart'] = 10000

        Jc = derivative(F_u_corrm, unc, u)
        ucproblem = NonlinearVariationalProblem(F_u_corrm, unc, bcs = list(bcu), J = Jc)
        ucsolver = NonlinearVariationalSolver(ucproblem)
        prmc = ucsolver.parameters
        prmc["newton_solver"]["absolute_tolerance"] = 1E-6
        prmc["newton_solver"]["relative_tolerance"] = 1E-6
        prmc["newton_solver"]["maximum_iterations"] = 50
        prmc["newton_solver"]["relaxation_parameter"] = .99
        prmc["newton_solver"]["error_on_nonconvergence"] = True
        prmc['linear_solver'] = 'gmres'
        prmc['preconditioner'] = 'hypre_euclid'
#        prmc['linear_solver'] = 'bicgstab'
#        prmc['preconditioner'] = 'jacobi'
        prmc['krylov_solver']['preconditioner']['reuse'] = True
        prmc['krylov_solver']['absolute_tolerance'] = 1E-7
        prmc['krylov_solver']['relative_tolerance'] = 1E-6
        prmc['krylov_solver']['maximum_iterations'] = 20000
        prmc['krylov_solver']['monitor_convergence'] = True
        prmc['krylov_solver']['nonzero_initial_guess'] = True
        prmc['krylov_solver']['gmres']['restart'] = 10000
        print 'Second Newton solver ok'


        # Time loop
        self.start_timing()
        while t<T:
	    t+=dt		
            self.timer("update & fetch bc")
            set_log_active(False)

            # Compute tentative velocity step
            usolver.solve()
            prm["newton_solver"]["maximum_iterations"] = 20
            prm["newton_solver"]["relaxation_parameter"] = .99
            prm["newton_solver"]["error_on_nonconvergence"] = True
            # Pressure correction I
            b = assemble(L_p_corr)
            if len(bcp) == 0:# or is_periodic(bcp):
                solver_p = solver_p_periodic
                normalize(b)
            else:
                solver_p = solver_p_dirichlet
            for bc in bcp: bc.apply(A_p_corr, b)
            solver_p_corr.solve(A_p_corr, p1.vector(), b)
            if len(bcp) == 0:# or is_periodic(bcp): 
		normalize(p1.vector())

            # Velocity correction I
            b = assemble(L_u_corr)
            for bc in bcu: bc.apply(A_u_corr, b)
            solver_u_corr.solve(A_u_corr, u2.vector(), b)


            # Pressure correction II
            b = assemble(L_p_corrII)
            if len(bcp) == 0:# or is_periodic(bcp):
                solver_p = solver_p_periodic
                normalize(b)
            else:
                solver_p = solver_p_dirichlet
            for bc in bcp: bc.apply(A_p_corrII, b)
            solver_p_corr.solve(A_p_corrII, p1.vector(), b)
            if len(bcp) == 0:# or is_periodic(bcp): 
		normalize(p1.vector())

            # Velocity correction momentum eq
            unc.assign(u2)
            ucsolver.solve()

            self.update(problem, t, unc, p1)
            u0.assign(unc)
            p0.assign(p1)
            prm['linear_solver'] = 'bicgstab'
            prm['preconditioner'] = 'jacobi'
            prmc['linear_solver'] = 'bicgstab'
            prmc['preconditioner'] = 'jacobi'
        return unc, p1

    def __str__(self):
        return "PISO"
