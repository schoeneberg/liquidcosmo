import numpy as np
from .folder import folder
from .matplotlib_defaults import default_settings
class foldercollection:
  
  def __init__(self):
    self.folderlist = []
    
  @classmethod
  def load(obj, *args, **kwargs):
    a = obj()
    if len(args) == 1:
      newargs = args[0]
    else:
      newargs = list(args)
    a.folderlist = [folder.load(arg,**kwargs) for arg in newargs]
    return a

  @property
  def folen(self):
    return len(self.folderlist)
  @property
  def bestfit(self):
    return [f.bestfit for f in self.folderlist]
  @property
  def N(self):
    return [f.N for f in self.folderlist]
  @property
  def chain(self):
    return [f.chain for f in self.folderlist]
  def copy(self):
    a = foldercollection()
    a.folderlist = [f.copy() for f in self.folderlist]
    return a
  def deepcopy(self):
    a = foldercollection()
    a.folderlist = [f.deepcopy() for f in self.folderlist]
    return a
  @property
  def cosmoargs(self):
    return [f.cosmoargs() for f in self.folderlist]
  @property
  def cosmopars(self):
    return [f.cosmopars() for f in self.folderlist]
  def cov(self):
    return [f.cov() for f in self.folderlist]
  def mean(self):
    return [f.mean() for f in self.folderlist]
  def names(self):
    return [f.names() for f in self.folderlist]
  def common_names(self):
    return list(set([f.names() for f in self.folderlist]))
  @property
  def d(self):
    return [f.d for f in self.folderlist]
  def derive(self, name, func, texname = None):
    flag = False
    for f in self.folderlist:
      try:
        f.derive(name, func, texname = texname)
        flag = True
      except KeyError as e:
        pass
    if not flag:
      raise Exception("Could not derive the asked parameter '{}' within any of the underlying folders. Make sure that the function/array you are passing as the 'func' parameter is correct for at least one of the underlying folders.".format(name))
  def get_chain(self,excludesmall=True,burnin_threshold=3):
    return [f.get_chain(excludesmall=excludesmall,burnin_threshold=burnin_threshold) for f in self.folderlist]
  def get_masked(self, mask):
    a = self.deepcopy()
    for i in range(a.folen):
      a.folderlist[i] = a.folderlist[i].get_masked(mask) 
    return a
  @property
  def logfile(self):
    return [f.logfile for f in self.folderlist]
  def set_range(self,parname,lower=None,upper=None):
    for f in self.folderlist:
      f.set_range(parname, lower=lower, upper=upper)
  def set_texname(self,parname,texname):
    flag = False
    for f in self.folderlist:
      try:
        f.set_texname(parname, texname)
        flag = True
      except:
        pass
    if not flag:
      raise Exception("Parameter '{}' not found in any of the folders contained in this collection.".format(parname))
  def subrange(self, parname=None, value=None):
    for f in self.folderlist:
      f.subrange(parname=parname, value=value)
    return self
  def write(self, fnames):
    assert(len(fnames) == self.folen)
    for i,f in enumerate(self.folderlist):
      f.write(fnames[i])
  @property
  def samples(self):
    return [f.samples for f in self.folderlist]
  @property
  def lens(self):
    return [f.lens for f in self.folderlist]
  def _readjust_bounds(self):
    res = self.copy()
    fbounds = [f.get_bounds() for f in res.folderlist]
    bounds = {}
    for bound in fbounds:
      for par in bound:
        if par in bounds:
          minbound = min((x for x in [bound[par][0],bounds[par][0]] if x is not None), default=None)
          maxbound = max((x for x in [bound[par][1],bounds[par][1]] if x is not None), default=None)
          bounds[par] = [minbound,maxbound]
        else:
          bounds[par] = bound[par]
    for f in res.folderlist:
      for par in bounds:
        f.set_range(par,lower=bounds[par][0],upper=bounds[par][1])
    return res
  def plot_getdist(self,ax=None,colors=None,alphas=None,**kwargs):
    from getdist.plots import get_subplot_plotter
    res = self._readjust_bounds()
    gdfolders = [f.to_getdist() for f in res.folderlist]
    spp = get_subplot_plotter(settings=default_settings)
    if 'filled' not in kwargs:
      kwargs['filled']=True
    spp.triangle_plot(gdfolders,alphas=alphas,colors=colors,
      line_args=([{'color':c} for c in colors] if colors else None),**kwargs)
    return spp
  def to_getdist(self):
    return [f.to_getdist() for f in self.folderlist]
  def __getitem__(self,q):
    if isinstance(q,(int,np.integer)):
      return [f[q] for f in self.folderlist]
    elif not isinstance(q,str):
      res = foldercollection()
      flag = False
      for f in self.folderlist:
        commonset = set(q).intersection(set(f.names))
        if commonset:
          res.folderlist.append(f[list(commonset)])
          flag = True
      if not flag:
        raise Exception("Could not find any of '{}' in any of the folders contained in this collection.".format(q))
      return res
    else:
      res = []
      flag = False
      for f in self.folderlist:
        if q not in f.names:
          continue
        else:
          res.append(f[q])
          flag = True
      if not flag:
        raise Exception("Could not find '{}' in any of the folders contained in this collection.".format(q))
      return np.array(res,dtype=object)
  def keep_only(self,*q):
    if isinstance(q,(int,np.integer)) or isinstance(q,str):
      return self[q]
    else:
      res = foldercollection()
      flag = False
      for f in self.folderlist:
        try:
          res.folderlist.append(f[q])
          flag = True
        except:
          pass
      if not flag:
        raise Exception("Could not find any of '{}' in any of the folders contained in this collection.".format(q))
      return res
  def __setitem__(self,q,v):
    if isinstance(v,foldercollection):
      for vi, vf in enumerate(v.folderlist):
        self.folderlist[vi][q] = vf
    elif (isinstance(v,list) or type(v) is np.ndarray) and len(v)>0 and (type(v[0]) is np.ndarray) and len(v)==self.folen:
      for i, f in enumerate(self.folderlist):
        f[q] = v[i]
    else:
      for f in self.folderlist:
        f[q] = v
  def __str__(self):
    return "Folderlist"+str(self.folderlist)

