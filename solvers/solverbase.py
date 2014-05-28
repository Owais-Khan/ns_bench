__author__ = "Anders Logg <logg@simula.no>"
__date__ = "2008-04-03"
__copyright__ = "Copyright (C) 2008-2010 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from dolfin import *
import sys
from time import time
import os
from os import getpid
from commands import getoutput
from numpy import linspace

# Common solver parameters
maxiter = default_maxiter = 200
tolerance = default_tolerance = 1e-4

master = MPI.process_number() == 0

class SolverBase:
    "Base class for all solvers."

    def __init__(self, options):

        # Store options
        self.options = options
	self.options["current_timestep"]=0

        # Reset some solver variables
        self._time = None
        self._cputime = 0.0
        self._timestep = 0


        #Keep track when to restart
        self.WCT=time()
        self.time_out=self.WCT+self.options["time_out"]*60.


        #If we want to restart the simulation
        if self.options["restart"]==True:
                 self._timestep = float(self.options["restart_time"])+1
	
	#Keep track of when peak systole occurs
	self.Qsys_counter=0

        # Reset files for storing solution
        self._ufiles = None
        self._pfile  = None

        # Reset storage for functional values and errors
        self._t = []
        self._M = []
        self._m = []
        self._e = []

        self._timer = None
        self._func_cache = {}

    def getMyMemoryUsage(self):
        mypid = getpid()
        mymemory = getoutput("ps -o rss %s" % mypid).split()[1]
        return mymemory

    def timer(self, msg):
        if self.options['timer'] and master:
            print "%10.0f ms: %s"%((time()-self._timer)*1000, msg)
            self._timer = time()

    def start_timing(self):
        """Start timing, will be paused automatically during update
        and stopped when the end-time is reached."""
        self._time = time()
        self._timer = time()

    def solve(self, problem, dt, plot_solution=True):
        "Solve problem"
        raise NotImplementedError

    def prefix(self, problem):
        "Return file prefix for output files"
        p = problem.__module__.split(".")[-1].lower()
        s = self.__module__.split(".")[-1].lower()
        return problem.output_location + p + "_" + s

    def desegregate(self, u):
        if not isinstance(u, (list, tuple)):
            return u
        if len(u) == 1:
            return u[0]
        assert MPI.num_processes() == 1
        V = u[0].function_space()
        W = VectorFunctionSpace(V.mesh(), V.ufl_element().family(), V.ufl_element().degree())
        f = self._func_cache.setdefault(str(W), Function(W))
        fv = f.vector()
        for i in range(len(u)):
            uv = u[i].vector()
            # FIXME: Assumptions about ordering of vector-function entries
            fv[i*len(uv):(i+1)*len(uv)] = uv
        return f

    def update(self, problem, t, u, p):
        "Update problem at time t"

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        # Update problem FIXME: Should this be called before problem.functional??
        problem.update_problem(t, u, p)

        # Evaluate functional and error
        m = problem.reference(t)
        M = problem.functional(t, u, p)
        if m is None:
            e = None
            print "M = %g (missing reference value)" % M
        else:
            e = abs(M - m)
            print "M = %g (reference %g), error = %g (maximum %g)" % (M, m, e, max([e] + self._e))

        # Store values
        self._t.append(t)
        self._M.append(M)
        self._m.append(m)
        self._e.append(e)
    def update_temp(self, problem, t, u, p):
        "Update problem at time t"
	self.options["current_timestep"]=self._timestep

	#Determine the current cycle Numer
        if self.options["current_cycle"]!=None or self.options["restart"]==True:
              if self._timestep%self.options["timesteps"]==0 and self._timestep!=0:
                      self.options["current_cycle"]+=1
		      self.options["file_name"]=self.options["case_name"]+"_curcyc_"+str(self.options["current_cycle"])+"_"
	if self.options["file_name"]==None:
		self.options["file_name"]=self.options["case_name"]

	filename=self.options["file_name"]+"t_"+str(round(t,2))+"_"

	#elif self.options["file_name"]==None:
	#	self.options["file_name"]=self.options["case_name"]
	#	filename=self.options["file_name"]
	#else:
	  #      self.options["file_name"]=self.options["case_name"]+"_curcyc_"+str(self.options["current_cycle"])+"_"
	   #     filename=self.options["file_name"]
			
	casedir = os.path.join("results",self.options["case_name"])
	casedir_check_point= os.path.join("results",self.options["case_name"],"check_point")
        u = as_object(u)

        if self.options['segregated']:
           s = max(ui.vector().norm('linf') for ui in u)
        else:
           s = u.vector().norm('linf')
        if s > 5 * getattr(problem, 'U', float('inf')):
            warning("A component in u is %.4g times characteristic velocity U"%round(s))
        if s > 1e10:
            raise RuntimeError("Runaway solution")

        # Add to accumulated CPU time
        timestep_cputime = time() - self._time
        self._cputime += timestep_cputime

        # Compute divergence
        if self.options["compute_divergence"]:
            check_divergence(u, p.function_space())

        # Update problem FIXME: Should this be called before problem.functional??
        problem.update_problem(t, u, p)

        # Evaluate functional and error
        m = problem.reference(t)
        M = problem.functional(t, u, p)
        if m is None:
             e = None
        #    if master:
        #        print "M = %s (missing reference value)" % fmt(M)
        else:
            e = abs(M - m)
            if master:
                print "M = %s (reference %s), error = %s (maximum %g)" % (fmt(M), fmt(m), fmt(e), max([e] + self._e))



        # Store values
                                                                                 
	self._t.append(t)
        self._M.append(M)
        self._m.append(m)
        self._e.append(e)




        # Save solution
        if self.options["save_solution"]:
            if (self._timestep - 1) % self.options["save_frequency"] == 0:
                # Create files
                if self._ufiles is None:
                    if self.options['segregated']:
                        # added '_' to separate i from vtu numbering
                        self._ufiles = [File(os.path.join(casedir, filename+"u%d_.pvd" % i), "compressed")
                                        for i in range(len(u))]
                    else:
                        self._ufiles = File(os.path.join(casedir, filename+"u.pvd"), "compressed")
                if self._pfile is None:
                    self._pfile = File(os.path.join(casedir, filename+"p.pvd"), "compressed")

                # Write to files
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        self._ufiles[i] << ui
                else:
                    self._ufiles << u
                self._pfile << p

        # Save solution at t = T
        if self.options["save_solution_at_t=T"]:
            if t >= problem.T:
                # Create files
                if self._ufiles is None:
                    self._ufiles = [File(os.path.join(casedir, filename+"u%d_at_end.pvd" % i), "compressed")
                                    for i in range(len(u))]
                else:
                    self._ufiles = File(os.path.join(casedir, filename+"u_at_end.pvd"), "compressed")
                if self._pfile is None:
                    self._pfile = File(os.path.join(casedir, filename+"p_at_end.pvd"), "compressed")

                # Write to files
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        self._ufiles[i] << ui
                else:
                    self._ufiles << u
                self._pfile << p



        # Save vectors in xml format
        if self.options["save_xml"]:
            if (self._timestep) % self.options["save_frequency"] == 0:
                timestr = "at_t%d_%.6e" % (self._timestep, t)
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        file = File(os.path.join(casedir, filename+"u%d_ts=%05.0f.xml.gz"% (i, self._timestep)))
                        file << ui #.vector()
                else:
                    file = File(os.path.join(casedir, filename+"u0_ts=%05.0f.xml.gz" % (self._timestep)))
                    file << u #.vector() #u[0].vector()
                file = File(os.path.join(casedir, filename+"p_ts=%05.0f.xml.gz" % (self._timestep,)))
                file << p


        # Save vector at check point
        if self.options["save_check_points"]:
            if (self._timestep - 1) % self.options["check_point_frequency"] == 0:
                timestr = "at_t%d_%.6e" % (self._timestep, t)
                if self.options['segregated']:
                    for i, ui in enumerate(u):
                        file = File(os.path.join(casedir_check_point, filename+"u%d_ts=%05.0f.xml.gz"% (i, self._timestep)))
                        file << ui #.vector()
                else:
                    file = File(os.path.join(casedir_check_point, filename+"u0_ts=%05.0f.xml.gz" % (self._timestep)))
                    file << u #.vector() #u[0].vector()
                file = File(os.path.join(casedir_check_point, filename+"p_ts=%05.0f.xml.gz" % (self._timestep,)))
                file << p
	   
	    #Save Last time step 
            if self.options["save_previous_ts"]==True:
                    if (self._timestep) % self.options["check_point_frequency"] == 0:
                        timestr = "at_t%d_%.6e" % (self._timestep, t)
                        if self.options['segregated']:
                            for i, ui in enumerate(u):
                                file = File(os.path.join(casedir_check_point, filename+"u%d_ts=%05.0f.xml.gz"% (i, self._timestep)))
                                file << ui #.vector()
                        else:
                            file = File(os.path.join(casedir_check_point, filename+"u0_ts=%05.0f.xml.gz" % (self._timestep)))
                            file << u #.vector() #u[0].vector()
                        file = File(os.path.join(casedir_check_point, filename+"p_ts=%05.0f.xml.gz" % (self._timestep,)))
                        file << p


        # Plot solution
        if self.options["plot_solution"]:
            # Plot velocity and pressure
            plot(self.desegregate(u), title="Velocity", rescale=True)
            plot(p, title="Pressure", rescale=True)

        # Check memory usage
        if self.options["check_mem_usage"]:
            if (self._timestep - 1) % self.options["check_frequency"] == 0:
                print 'Memory usage is:' , self.getMyMemoryUsage()

        # Print progress
        if master:
            ss = self._cputime * (problem.T/t-1)
            hh, ss = divmod(ss, 60*60)
            mm, ss = divmod(ss, 60)
            print
            s = "Time step %d finished in %.2f seconds, %.1f%% done (t=%.3g, T=%g; %02d:%02d:%02d remaining)." \
                % (self._timestep, timestep_cputime, 100.0*(t / problem.T), t, problem.T, hh, mm, ss)
            print "-"*len(s)
	    print s
            print "-"*len(s)


        # Increase time step and record current time
        self._timestep += 1
        self._time = time()




        ####AUTOMATED SIMULATION RESTART###     
        #If simulation is about to end, restart the simulation
        if time()>=self.time_out:
                #Make restart option true
                self.options["restart"]=True

                #List the files in directory according to their update time 
                #and obtain the latest file
                restart_path="./results/"+self.options["case_name"]+"/check_point/"
                files=sorted(os.listdir(restart_path),key=lambda p: os.path.getctime(os.path.join(restart_path,p)))



                #check to see which files are the velocity/pressure files
                files_new=[]
                for i in files:
                        test_file=i.find("ts=")
                        if test_file>0:
                                files_new.append(i)

                restart_time=int(files_new[-1].split("=")[1].split(".xml.gz")[0])


                if restart_time%self.options["save_frequency"]==0:
                        restart_time=restart_time-(self.options["save_frequency"]-1)


                self.options["restart_time"]=restart_time
                #self.options["restart_time"]=int(files[-1].split("=")[1].split(".xml.gz")[0])

                #update the no of restart
                self.options["no_of_restart"]+=1

                #Print message
                if master:
                        cmd=self.run_me()

                        print "\n Maximum Wall-Clock-Time reached"
                        print "Restarting simulation at Time Step = ", self.options["restart_time"]
                        print "executing: ",cmd

                        #Restart the simulation with new data   
                        cmd=self.run_me()
                        os.system(cmd)
                        sys.exit()


    def run_me(self):
        #First Piece
        command1="""#!/bin/bash 
# MOAB/Torque submission script for SciNet GPC
#
#PBS -l nodes=%d:ppn=8,walltime=%02.0f:%02.0f:00"""%(self.options["processors"]/8,self.options["wct_hrs"],self.options["wct_mins"])
	
	command2 = "\n#PBS -N "+self.options["case_name"]+"\n"

        #Second Piece
        if self.options["node_type"]=="debug":
                command3="\n#PBS -q debug\n"
        else:
                command3="\n# #PBS -q debug\n"

        #Third Piece
        command4="""# DIRECTORY TO RUN - $PBS_O_WORKDIR is directory job was submitted from
ssh gpc03
cd $PBS_O_WORKDIR
#export INSTANT_CACHE_DIR=Work/instant-cache/
#export INSTANT_ERROR_DIR=Work/instant-error/
# EXECUTION COMMAND; -np = nodes*ppn
# mpirun --tag-output -np 8 env  | grep LD_LIBRARY_PATH\n"""



        temp=""
        for i in self.options.keys():
                temp=temp+i+"="+str(self.options[i]).replace(" ","")+" "

        command5="mpirun -n "+str(self.options["processors"])+" python ns "+str(self.options["problem_name"])+ " "+str(self.options["solver_name"])+" "+temp+" &> logs/"+str(self.options["case_name"])+"_restart_"+str(self.options["no_of_restart"])



        #create a run file
        run_me_file  =self.options["case_name"]+".sh"
        #Outfile
        outfile=open(run_me_file,"w")
        outfile.write(command1+command2+command3+command4+command5)
        outfile.close()
        #Outscript
        pwd=os.getcwd()
        cmd="ssh gpc03 'cd "+pwd+" ; qsub ./"+run_me_file+";'"
        return cmd







    def select_timestep(self, problem):
        "Return time step and number of time steps for problem."

        # FIXME: This looks very complex, should be cleaned up

        T  = problem.T
        U  = getattr(problem, "U", float("nan"))
        nu = problem.nu
        h  = MPI.min(problem.mesh.hmin())

        try:
            # FIXME: This sequence of ifs make no sense. Clean up...

            if self.options["dt"]:
                if master: print "Using user supplied dt"
                dt = self.options["dt"]
                n = int(T / dt + 0.5)
                dt = T / n

            elif self.options["dt_division"] != 0 and getattr(problem, "dt", 0) > 0:
                if master: print 'Using problem.dt and time step refinements'
                dt = problem.dt / int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n

            # Use time step specified in problem if available
            elif getattr(problem, "dt", 0) > 0:
                if master: print 'Using problem.dt'
                dt = problem.dt
                n  = int(T / dt)

            # Otherwise, base time step on mesh size
            elif self.options["dt_division"] != 0:
                if master: print 'Computing time step according to stability criteria and time step refinements'
                dt = 0.25*h**2 / (U*(nu + h*U))
                dt /= int(sqrt(2)**self.options["dt_division"])
                n  = int(T / dt + 1.0)
                dt = T / n

            # Otherwise, base time step on mesh size
            else:
                if master: print 'Computing time step according to stability criteria'
                # dt =  0.25*h**2 / (U*(nu + h*U))
                dt =  0.2*(h / U)
                n  = int(T / dt + 1.0)
                dt = T / n

        finally:
            if dt != dt:
                warning("Unable to compute dt. Please set U in problem, or specify dt in options")

        # Compute range
        t_range = linspace(0,T,n+1)[1:] # FIXME: Comment out [1:] to run g2ref g2ref

        # Compute save_frequency if save_number is given
        if self.options["save_number"] is not None:
            frequency = int(1 + len(t_range) / (1+self.options["save_number"]))
            if master:
                print "Selected save frequency of %d (total %d snapshots)" % (frequency, len(t_range)//frequency)
            self.options["save_frequency"] = frequency

        # Report time step
        if master:
            print " "
            print 'Number of timesteps:' , len(t_range)
            print 'Size of timestep:' , dt
            print " "

        self._num_steps = len(t_range)
        return dt, t_range[0], t_range

    def eval(self):
        "Return last functional value and maximum error in functional value on [0, T]"

        # Plot values
        if self.options["plot_functional"]:
            from pylab import plot, xlabel, ylabel, grid, show
            plot(self._t, self._M)
            xlabel('t')
            ylabel('Functional')
            grid(True)
            show()

        # Return value
        if self._e[0] is None:
            return self._M[-1], None
        else:
            return self._M[-1], max([0.0] + self._e)

    def _list_or_function(self, u):
        "Return a single function if possible, else a list."
        if self.options['segregated'] or not isinstance(u, (list, tuple)):
            return u
        else:
            return u[0]


    def cputime(self):
        "Return accumulated CPU time."
        return self._cputime

def fmt(M, format="%g"):
    if hasattr(M, "__iter__"):
        return "[%s]"%", ".join(format%m for m in M)
    else:
        return format%M

def as_object(u):
    "Return a single object if possible, else a list."
    if not isinstance(u, (list, tuple)) or len(u) > 1:
        return u
    else:
        return u[0]

def epsilon(u):
    "Return symmetric gradient."
    return 0.5*(grad(u) + grad(u).T)

def sigma(u, p, nu):
    "Return stress tensor."
    return 2*nu*epsilon(u) - p*Identity(u.cell().d)

#def is_periodic(bcs):
 #   "Check if boundary conditions are periodic."
  #  return all(isinstance(bc, PeriodicBC) for bc in bcs)

def has_converged(r, iter, method, maxiter=default_maxiter, tolerance=default_tolerance):
    "Check if solution has converged."
    if master: print "Residual = %.3g"%r
    if r < tolerance:
        if master: print "%s iteration converged in %d iteration(s)." % (method, iter + 1)
        return True
    elif iter == maxiter - 1:
        raise RuntimeError, "%s iteration did not converge." % method
    return False

def check_divergence(u, Q):
    "Check divergence of velocity."

    # Compute L2 norm of divergence
    if master: print "||div u||_L2 =", norm(u, "Hdiv0")

    # Compute projection of div u into Q_0
    pdivu = project(div(u), Q)
    zero = Constant(Q.mesh(), 0.0)
    bc = DirichletBC(Q, zero, DomainBoundary())
    bc.apply(pdivu.vector())

    # Compute "weak" L2 norm of divergence
    if master: print "||div u||_w  =", sqrt(abs(assemble(pdivu*div(u)*dx, mesh=Q.mesh())))
