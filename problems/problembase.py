
__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-03-19"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

# Modified by Kent-Andre Mardal, 2008.

from dolfin import *
from numpy import linspace
from math import *

master = MPI.process_number() == 0

class ProblemBase:
    "Base class for all problems."

    def __init__(self, options):

        # Store options
        self.options = options

        # Parameters must be defined by subclass
        self.mesh = None
        self.f    = None
        self.bcu  = []
        self.bcp  = []
        self.nu   = None
        self.t    = 0
        self.T    = None
        self.dt   = None
        self.u0   = None
        self.p0   = None
        self.u    = None
        self.p    = None
        self.U    = 1.0
        self.output_location = ''

        # FIXME: Temporary while testing
        self.bcpsi = []

    def update_problem(self, t, u, p):
        "Update problem at time t"

        # Update state
        self.t = t
        self.u = u
        self.p = p

        # Call problem-specific update
        self.update(t, u, p)

    def update(self, t, u, p):
        "Problem-speficic update at time t"
        pass

    def functional(self, t, u, p):
        "Return value of functional of interest"
        return 0.0

    def reference(self, t):
        "Return reference value for functional"
        return None

    def tolerance(self, problem):
        "Return tolerance (used as local convergence criterion)."
        if str(problem) == 'Channel':
            return 1e-11
        elif str(problem) == 'Cylinder':
            return 1e-7
        else:
            return 1e-6

    def uConstant(self, values):
        if isinstance(values, tuple) and self.options['segregated']:
            return [Constant(v) for v in values]

    def uExpr(self, cppcode, **kwargs):
        if self.options['segregated'] and isinstance(cppcode, (tuple, list)):
            return [Expression(e, **kwargs) for e in cppcode]
        else:
            return [Expression(cppcode, **kwargs)]



    def eval(self, func, point, gather=True):
        """Parallel-safe function evaluation"""
        if gather:
            if hasattr(func, 'update'):
                func.update() # dolfin dev
            else:
                func.gather() # dolfin 1.0
        if len(func.shape())==1:
            M = [0]*func.shape()[0]
        else:
            M = 0
        try:
            M = func(point)
            N = MPI.sum(1) # Succeeding processors participate in the MPI collective here
        except RuntimeError:
            N = MPI.sum(0) # Failing processors participate in the MPI collective here
            if N == 0:
                raise      # All processors failed
        if master and N > 1:
            warning("%d processors returned function value, which is unexpected (but probably ok)"%N)
        if hasattr(M, '__iter__'):
            for i in range(len(M)):
                M[i] = MPI.sum(M[i])/N
        else:
            M = MPI.sum(M)/N
        return M


    def uEval(self, func, component, point):
        if self.options['segregated']:
            return self.eval(func[component], point)
        else:
            return self.eval(func, point)[component]


    def timestep(self, problem):
        "Return time step and number of time steps for problem."

        # FIXME: This looks very complex, should be cleaned up

        T = problem.T
        U = problem.U
        nu =problem.nu
        h  = self.mesh.hmin()
        "Return time step and number of time steps for problem. Used for debugging / compilation only"
        if self.options["max_steps"] is not None:
            dt =  0.25*h**2 / (U*(nu + h*U))
            n  = self.options["max_steps"]
            T  = n*dt
            t_range = linspace(0, T, n + 1)[1:]
        else:

            if self.options["dt_division"] != 0 and problem.dt is not None:
                dt =  (problem.dt)/int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n
                print 'Using problem.dt and time step refinements'

            # Use time step specified in problem if available
            elif not problem.dt is None:
                dt = problem.dt
                print dt
                n  = int(T / dt)
                print 'Using problem.dt'

            # Otherwise, base time step on mesh size
            elif self.options["dt_division"] != 0:
                dt =  0.25*h**2 / (U*(nu + h*U))/int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n
                print 'Computing time step according to stability criteria and time step refinements'

            # Otherwise, base time step on mesh size
            else:
#                dt =  0.25*h**2 / (U*(nu + h*U))
                dt =  0.2*(h / U)
                n  = int(T / dt + 1.0)
                dt = T / n
                print 'Computing time step according to stability criteria'

        # Compute range
        t_range = linspace(0,T,n+1)[1:] # FIXME: Comment out [1:] to run g2ref g2ref

        # Report time step
        print " "
        print 'Number of timesteps:' , len(t_range)
        print 'Size of timestep:' , dt
        print " "

        return dt, t_range[0], t_range

    # FIXME: Remove this and merge with boundary_condition()
    def pressure_bc(self, Q):
        return Constant(0)

