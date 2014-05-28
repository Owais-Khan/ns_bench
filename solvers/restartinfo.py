__author__ = "Joachim B Haga <jobh@simula.no>"
__date__ = "2012-02-17"
__copyright__ = "Copyright (C) 2012 " + __author__
__license__  = "GNU GPL version 3 or any later version"

from dolfin import *
from numpy import linspace
import os
from glob import glob
import re

master = MPI.process_number() == 0

class RestartInfo(object):
    def __init__(self, casename, timestep=None, t=None):
        self.path = os.path.join("results", casename)
        regexp = re.compile("_at_t([0-9]*)_(.*).xml.gz")
        if timestep is not None:
            filenames = glob(os.path.join(self.path, "*_at_t%d_*.xml.gz"%timestep))
            m = filenames and regexp.search(filenames[0])
            if m:
                self.timestep = timestep
                self.t = float(m.group(2))
        elif t is not None:
            filenames = glob(os.path.join(self.path, "*_at_t*.xml.gz"))
            for fn in filenames:
                m = regexp.search(fn)
                if m and abs(float(m.group(2))-t) < abs(1e-5*t):
                    self.t = t
                    self.timestep = int(m.group(1))
                    break
        if not hasattr(self, 't'):
            raise RuntimeError("Unable to match time/timestep (no saved data in casename %s)"%casename)

    def select_timestep(self, dt, T):
        """Return new time step info"""
        n = int((T-self.t) / dt + 0.5)
        if n < 1:
            raise RuntimeError("Trying to restart from t >= T")
        self.t_range = linspace(self.t, T, n + 1)[1:]
        return (T-self.t)/n, self.t, self.t_range

    def p(self, t, Q):
        return self._load_at_time(t, Q, "p")

    def u(self, t, V):
        if V.num_sub_spaces() > 0:
            # Vector function space, not segregated
            return [self._load_at_time(t, V, "u")]
        else:
            # Scalar function space, segregated
            return [self._load_at_time(t, V, "u%d"%d)
                    for d in range(V.mesh().topology().dim())]

    def _load_at_time(self, t, functionspace, name):
        filenames = glob(os.path.join(self.path, "%s_at_t%d_*.xml.gz"%(name, self.timestep)))
        if len(filenames) != 1:
            if master:
                print filenames
            raise RuntimeError()
        regexp = re.compile("_at_t%d_(.*).xml.gz"%self.timestep)

        # FIXME: Should support going back in time, using linear interpolation
        # between saved files.
        if master and t != self.t:
            warning("Using data for current time %.2g instead of %.2g for %s"%(self.t, t, name))

        f = Function(functionspace)
        File(filenames[0]) >> f.vector()
        return f
