import numpy as np
import numbers
_round_exponential_max = 5
_round_exponential_min = -5
def round_reasonable(val, errp=None,errm=None , digits=1, equalize=0.0):
  # Do some accounting a priori
  if digits<0:
    raise ValueError("Cannot print negative number of significant digits")
  if errm is not None and errp is None:
    raise ValueError("Cannot give only negative but no positive error bar")
  # Check for arrays instead of values, convert individually
  if isinstance(val,(list,tuple,np.ndarray)):
    if errp and not isinstance(errp,(list,tuple,np.ndarray)) and not len(errp)==len(val):
      raise ValueError("Cannot pass a list of vals but only a single errp")
    if errm and not isinstance(errm,(list,tuple,np.ndarray)) and not len(errp)==len(val):
      raise ValueError("Cannot pass a list of vals but only a single errm")
    return [round_reasonable(v,errp=(errp[iv] if errp else None),errm=(errm[iv] if errm else None),digits=digits) for iv,v in enumerate(val)]
  # Only value, no errors
  if errp is None:
    if val==0:
      first_sigfig = 0
    else:
      first_sigfig = int(np.floor(np.log10(np.abs(val))))
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      return "{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \cdot 10^{"+str(first_sigfig)+"}"
    else:
      return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)
  # Symmetric error
  elif errm is None:
    if errp<0:
      raise ValueError("Cannot pass negative errp")
    # error = 0 is equivalent to no error
    if errp==0:
      return round_reasonable(val, errp=None,errm=None , digits=digits)
    # Use error sigfig to do rest
    first_sigfig = int(np.floor(np.log10(np.abs(errp))))
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      return "("+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \pm "+"{:.{acc}f}".format(errp*10**(-first_sigfig),acc=digits)+") \cdot 10^{"+str(first_sigfig)+"}"
    else:
      return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)+" \pm "+"{:.{acc}f}".format(errp,acc=-first_sigfig+digits)
  else:
    if errp<0:
      raise ValueError("Cannot pass negative errp")
    if errm<0:
      errm=-errm
    if errp==0 and errm==0:
      return round_reasonable(val, errp=None,errm=None , digits=digits)
    if errp==0:
      errp_sigfig = 0
    else:
      errp_sigfig = int(np.floor(np.log10(np.abs(errp))))
    if errm==0:
      errm_sigfig = 0
    else:
      errm_sigfig = int(np.floor(np.log10(np.abs(errm))))
    # Choose the LARGER error for the significant digit
    first_sigfig = max(errp_sigfig,errm_sigfig)
    if equalize != 0.:
      if np.abs(errm*10**(-first_sigfig+digits)-errp*10**(-first_sigfig+digits))<equalize:
        return round_reasonable(val, errp=(errp+errm)/2.,errm=(errp+errm)/2. , digits=digits, equalize=0.)
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      pluserrstr = "{:.{acc}f}".format(errp*10**(-first_sigfig),acc=digits)
      minuserrstr = "{:.{acc}f}".format(errm*10**(-first_sigfig),acc=digits)
      if pluserrstr!=minuserrstr:
        return "({"+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+"}^{+"+pluserrstr+"}_{-"+minuserrstr+"}"+") \cdot 10^{"+str(first_sigfig)+"}"
      else:
        return "("+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \pm "+pluserrstr+") \cdot 10^{"+str(first_sigfig)+"}"
    else:
      pluserrstr = "{:.{acc}f}".format(errp,acc=-first_sigfig+digits)
      minuserrstr = "{:.{acc}f}".format(errm,acc=-first_sigfig+digits)
      if pluserrstr!=minuserrstr:
        return "{"+"{:.{acc}f}".format(val,acc=-first_sigfig+digits)+"}^{+"+pluserrstr+"}_{-"+minuserrstr+"}"
      else:
        return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)+" \pm "+"{:.{acc}f}".format(errp,acc=-first_sigfig+digits)

# Since by default lists of arrays don't permit many operations, we define them here
class RaggedArray:
  def __init__(self, data):
    self.data = np.fromiter((np.array(d) for d in data),dtype=object)

  def __checked_op__(self, a, opname):
    if isinstance(a, RaggedArray):
      self.__check_compatible(a)
      return RaggedArray(np.fromiter((getattr(d, opname)(a.data[i]) for i,d in enumerate(self.data)),dtype=object))
    elif isinstance(a, numbers.Number):
      return RaggedArray(np.fromiter((getattr(d, opname)(a) for d in self.data),dtype=object))
    else:
      raise ValueError("Unknown operation '{}' for RaggedArray and {}".format(opname,a.__class__()))

  def __abs__(self):
    return RaggedArray(np.fromiter((d.__abs__() for d in self.data),dtype=object))
  def __add__(self, a):
    return self.__checked_op__(a, '__add__')
  def __ceil__(self):
    return RaggedArray(np.fromiter((d.__ceil__() for d in self.data),dtype=object))
  def __floor__(self):
    return RaggedArray(np.fromiter((d.__floor__() for d in self.data),dtype=object))
  def __floordiv__(self, a):
    return self.__checked_op__(a, '__floordiv__')
  def __divmod__(self, a):
    return self.__checked_op__(a, '__divmod__')
  def __eq__(self, a):
    return self.__checked_op__(a, '__eq__')
  def __ge__(self, a):
    return self.__checked_op__(a, '__ge__')
  def __gt__(self, a):
    return self.__checked_op__(a, '__gt__')
  def __le__(self, a):
    return self.__checked_op__(a, '__le__')
  def __lt__(self, a):
    return self.__checked_op__(a, '__lt__')
  def __mod__(self, a):
    return self.__checked_op__(a, '__mod__')
  def __mul__(self, a):
    return self.__checked_op__(a, '__mul__')
  def __ne__(self, a):
    return self.__checked_op__(a, '__ne__')
  def __neq__(self):
    return RaggedArray(np.fromiter((d.__neg__() for d in self.data),dtype=object))
  def __pow__(self, a):
    return self.__checked_op__(a, '__pow__')
  def __rsub__(self, a):
    return self.__checked_op__(a, '__rsub__')
  def __rmul__(self, a):
    return self.__checked_op__(a, '__rmul__')
  def __rdiv__(self, a):
    return self.__checked_op__(a, '__rdiv__')
  def __rtruediv__(self, a):
    return self.__checked_op__(a, '__rtruediv__')
  def __rfloordiv__(self, a):
    return self.__checked_op__(a, '__rfloordiv__')
  def __rmod__(self, a):
    return self.__checked_op__(a, '__rmod__')
  def __rdivmod__(self, a):
    return self.__checked_op__(a, '__rdivmod__')
  def __rpow__(self, a):
    return self.__checked_op__(a, '__rpow__')
  def __rlshift__(self, a):
    return self.__checked_op__(a, '__rlshift__')
  def __rrshift__(self, a):
    return self.__checked_op__(a, '__rrshift__')
  def __truediv__(self, a):
    return self.__checked_op__(a, '__truediv__')
  def __sub__(self, a):
    return self.__checked_op__(a, '__sub__')
  def __trunc__(self, a):
    return self.__checked_op__(a, '__trunc__')

  def __getitem__(self,q):
    return self.data[q]
  def __setitem__(self,q,v):
    self.data[q] = v
  def __repr__(self):
    return "RaggedArray[{}]".format(",".join([repr(d) for d in self.data]))
  def __str__(self):
    return "RaggedArray[{}]".format(",".join([str(d) for d in self.data]))
  def __check_compatible(self, a):
    if not isinstance(a,RaggedArray):
      raise Exception("Comparing compatibility with non-RaggedArray")
    if len(self.data)!=len(a.data):
      raise ValueError("Incompatible length of RaggedArray ({} vs {})".format(len(self),len(a)))
    broadcastable = True
    for i in range(len(self.data)):
      if not all((m == n) or (m == 1) or (n == 1) for m, n in zip(self.data[i].shape[::-1], a.data[i].shape[::-1])):
        raise ValueError("Dimension {} of RaggedArray is not compatible ({} vs {})".format(i, self.data[i].shape, a.data[i].shape))
