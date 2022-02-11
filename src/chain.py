from copy import deepcopy
import numpy as np
# Basically a wrapper for a dictionary of arrays, but with some added basic functionality
import matplotlib.pyplot as plt
import getdist
from collections import OrderedDict

class chain:

  # Initialize the chain object
  def __init__(self, dictionary):
    if not isinstance(dictionary,OrderedDict):
      raise Exception("Must supply an _ordered_ python dictionary to create a chain object")
    self._d = dictionary
    self.names = list(self._d.keys())
    self.N = len(self._d[self.names[0]])
    self._idx_bestfit = None

  # Get an item from the chain
  def __getitem__(self, q):
    if isinstance(q,(int,np.integer)):
      return chain(OrderedDict({key:[self._d[key][q]] for key in self.names}))
    elif isinstance(q,str):
      return self._d[q]
    elif isinstance(q,slice):
      return chain(OrderedDict({key:self._d[key][q] for key in self.names}))
    elif isinstance(q,np.ndarray) and q.dtype=="bool":
      return chain(OrderedDict({key:self._d[key][q] for key in self.names}))
    elif isinstance(q,list) and len(q)>0 and isinstance(q[0],str):
      return chain(OrderedDict({key:self._d[key] for key in q}))
    else:
      raise Exception("Cannot get from chain with object of type "+str(type(q)))
  def get_dict(self,i):
    return {key:self._d[key][i] for key in self.names}

  # Set an item in the chain
  def __setitem__(self,q,v):
    if isinstance(q,str):
      if (isinstance(v,(list,np.ndarray)) and len(v) == self.N):
        if q not in self.names:
          self.names.append(q)
        self._d[q] = v
      else:
        raise Exception("Cannot set part '{}' of chain with '{}' -- invalid type/dimension.".format(q,v))
    elif isinstance(q,slice):
      if isinstance(v,chain):
        for key in self.names:
          if key in v.names:
            self._d[key][q] = v[key]
        return chain(OrderedDict({key:self._d[key][q] for key in self.names}))
      else:
        raise Exception("Cannot set part of chain with a non-chain object of type "+str(type(v)))
    elif isinstance(q,int):
      return chain(OrderedDict({key:self._d[key][q] for key in self.names}))
    else:
      raise Exception("Cannot set chain with object of type "+str(type(q)))
  def __str__(self):
    return "Chain"+self._str_part()
  def _str_part(self):
    if self.N == 1:
      return "["+",".join(str(key) for key in self.names)+"|len="+str(self.N)+"] = ["+",".join(str(self._d[key][0]) for key in self.names)+"]"
    return "["+",".join(str(key) for key in self.names)+"|len="+str(self.N)+"]"

  @property
  def bestfit(self):
    if (self._idx_bestfit is None):
      self._idx_bestfit = np.argmin(self['lnp'])
    return self[self._idx_bestfit]
  def derive(self, name, func):
    self[name] = func(self)
      
