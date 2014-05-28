from dolfin import *
from ufl.form import Form

class RhsGenerator(object):
    """Class for storing the instructions to create the RHS vector b.
    The two main purposes of this class are:

    - make it easy to define the LHS and RHS in the same place

    - make it easy to generate RHS from matrix-XXX products, where XXX may be either
      * a Constant (which can be projected to a vector at once)
      * an Expression (which must be projected each time, because its parameters may change)
      * a Function
      """

    def __init__(self, space):
        self.space = space
        self.matvecs = []
        self.form = None
        self.vecs = []

    def __iadd__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            assert isinstance(A, GenericMatrix)
            self.matvecs.append((A, self._as_vector_or_timedep(x), 1))
        elif isinstance(ins, GenericVector):
            self.vecs.append(ins)
        elif isinstance(ins, Form):
            if self.form is None:
                self.form = ins
            else:
                self.form += ins
        else:
            raise RuntimeError, "Unknown RHS generator "+str(type(ins))
        return self

    def __isub__(self, ins):
        if isinstance(ins, tuple):
            A, x = ins
            if isinstance(A, GenericMatrix):
                self.matvecs.append((A, self._as_vector_or_timedep(x), -1))
                return self
        raise RuntimeError, "Try '+=' instead"

    def _as_vector_or_timedep(self, x):
        if isinstance(x, (GenericVector, Expression, Function)):
            return x
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def _as_vector(self, x):
        if isinstance(x, GenericVector):
            return x
        if isinstance(x, Function):
            return x.vector()
        return assemble(inner(x, TestFunction(self.space)) * dx)

    def __call__(self, bcs=None, symmetric_mod=None):
        f = Function(self.space)
        b = f.vector().copy() # dolfin bug 889021
        for mat, x, alpha in self.matvecs:
            b_ = mat * self._as_vector(x)
            if alpha != 1:
                b_ *= alpha
            b += b_
        for vec in self.vecs:
            b += vec
        if self.form:
            assemble(self.form, tensor=b, add_values=True, reset_sparsity=False)
        for bc in self._wrap_in_list(bcs, "bcs", DirichletBC):
            bc.apply(b)
        if symmetric_mod:
            b -= symmetric_mod*b
        return b

    def _wrap_in_list(self, obj, name, types=type):
        if obj is None:
            lst = []
        elif hasattr(obj, '__iter__'):
            lst = list(obj)
        else:
            lst = [obj]
        for obj in lst:
            if not isinstance(obj, types):
                raise TypeError("expected a (list of) %s as '%s' argument" % (str(types),name))
        return lst
