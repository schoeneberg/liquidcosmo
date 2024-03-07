import numpy as np
from .folder import folder
from .matplotlib_defaults import default_settings
from itertools import cycle
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
  @property
  def tags(self):
    return [f.tag for f in self.folderlist]
  def cov(self,parnames=None):
    return [f.cov(parnames=parnames) for f in self.folderlist]
  def mean(self,parnames=None,asdict=False):
    return [f.mean(parnames=parnames,asdict=asdict) for f in self.folderlist]
  def std(self,parnames=None,asdict=False):
    return [f.std(parnames=parnames,asdict=asdict) for f in self.folderlist]
  def credible(self,parnames=None,p=None,sigma=None,twoside=False,upper=True):
    return [f.credible(parnames=parnames,p=p,sigma=sigma,twoside=twoside,upper=upper) for f in self.folderlist]
  def constraint(self,parnames=None):
    return [f.constraint(parnames=parnames) for f in self.folderlist]
  def texconstraint(self,parnames=None,withname=True):
    return [f.texconstraint(parnames=parnames,withname=withname) for f in self.folderlist]
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
  def set_range(self,parname,lower=None,upper=None,destructive=False):
    for f in self.folderlist:
      f.set_range(parname, lower=lower, upper=upper,destructive=destructive)
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
  def gelman(self, parnames=None,subdivisions=None):
    return [f.gelman(parnames=parnames,subdivisions=subdivisions) for f in self.folderlist]
  def max_gelman(self, subdivisions=None):
    return [f.max_gelman(subdivisions=subdivisions) for f in self.folderlist]
  @property
  def samples(self):
    return [f.samples for f in self.folderlist]
  @property
  def lens(self):
    return [f.lens for f in self.folderlist]
  def merge(self,basechain=0):
    lensums = np.concatenate([[0],np.cumsum(self.lens)])
    tot_len = lensums[-1]
    obj = self.folderlist[basechain].deepcopy()
    obj.chain.N = tot_len
    for i,name in enumerate(obj.names):
      obj.chain[name] = np.empty(tot_len) # destroy it all first (keeping only metainfo)
    for j, fo in enumerate(self.folderlist):
      for i,name in enumerate(obj.names):
        obj.chain[name][lensums[j]:lensums[j+1]] = fo.chain[name]
    return obj
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
  def plot_getdist(self,ax=None,colors=None,alphas=None,add_point=None,contours=2,**kwargs):
    from getdist.plots import get_subplot_plotter
    res = self._readjust_bounds()
    gdfolders = [f.to_getdist() for f in res.folderlist]
    ana_set = kwargs.pop('analysis_settings',None)
    if ana_set is not None:
      for gdf in gdfolders:
        gdf.updateSettings(settings=ana_set)
    if 'filled' not in kwargs:
      kwargs['filled'] = True
    if 'legend_labels' not in kwargs:
      kwargs['legend_labels'] = [f.tag for f in res.folderlist]
    line_args = kwargs.pop('line_args',[{} for i in range(len(gdfolders))])
    if "linestyle" in kwargs:
      if isinstance(kwargs['linestyle'],(list,tuple,np.ndarray)):
        for i in range(len(gdfolders)):
          line_args[i].update({"ls":kwargs['linestyle'][i%len(kwargs['linestyle'])]})
      else:
        line_args = [kwargs['linestyle'] for i in range(len(gdfolders))]
    if not colors:
      cyc = cycle(default_settings.solid_colors)
      colors = [next(cyc) for i in range(len(gdfolders))]
    for i in range(len(gdfolders)):
      line_args[i].update({"color":colors[i%len(colors)]})
    contour_ls = kwargs.pop('contour_ls',[line_args[i].get('ls','-') for i in range(len(gdfolders))])
    spp = get_subplot_plotter(settings=default_settings,width_inch=kwargs.pop('width_inch',None),subplot_size_ratio=kwargs.pop('subplot_size_ratio',None))
    spp.settings.num_plot_contours = contours
    rect = kwargs.pop('rectangle',None)
    if rect != None:
      spp.rectangle_plot(rect['x'],rect['y'],roots=gdfolders, alphas=alphas,colors=colors,contour_ls=contour_ls,line_args=line_args,**kwargs)
    else:
      spp.triangle_plot(gdfolders, alphas=alphas,colors=colors,contour_ls=contour_ls,line_args=line_args,**kwargs)
    # Delegate to first folder, to use same function
    self.folderlist[0]._add_point(spp,add_point,names=self.common_names)
    return spp

  def to_getdist(self):
    return [f.to_getdist() for f in self.folderlist]
  def __getitem__(self,q):
    wants_sublist = (
        (isinstance(q,(tuple,list)) and all([isinstance(t,(int,np.integer)) for t in q]))
        or (isinstance(q,np.ndarray) and np.issubdtype(q.dtype, np.integer))
      )
    if isinstance(q,(int,np.integer)):
      res = foldercollection()
      res.folderlist = [self.folderlist[q]]
      return res
    elif isinstance(q,slice):
      res = foldercollection()
      res.folderlist = self.folderlist[q]
      return res
    elif wants_sublist:
      res = foldercollection()
      for i in q:
        if i>=len(self.folderlist) or i<0:
          raise IndexError("{} folders contained in this collection, but asked for element {}".format(len(self.folderlist),i))
      res.folderlist = [self.folderlist[i] for i in q]
      return res
    elif not isinstance(q,str):
      res = foldercollection()
      flag = False
      for f in self.folderlist:
        commonset = set(q).intersection(set(f.names))
        if commonset:
          res.folderlist.append(f[list(name for name in q if name in commonset)])
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

