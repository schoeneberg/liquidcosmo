from .folder import folder
class foldercollection:
  
  def __init__(self):
    self.folderlist = None
    
  @classmethod
  def load(obj, *args):
    a = obj()
    if len(args) == 1:
      newargs = args[0]
    else:
      newargs = list(args)
    a.folderlist = [folder.load(arg) for arg in args]
    return a

  @property
  def bestfit(self):
    return [f.bestfit for f in self.folderlist]
  @property
  def N(self):
    return [f.N for f in self.folderlist]
  @property
  def chain(self):
    return [f.chain for f in self.folderlist]
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
  def derive(self, name, func, verbose = 0, texname = None):
    for f in self.folderlist:
      f.derive(name, func, verbose = verbose, texname = texname)
  def get_chain(self,excludesmall=True,burnin_threshold=5):
    return [f.get_chain(excludesmall=excludesmall,burnin_threshold=burnin_threshold) for f in self.folderlist]
  def get_masked(self, mask):
    a = self.deepcopy()
    for i in range(len(a.folderlist)):
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
    assert(len(fnames) == len(self.folderlist))
    for i,f in enumerate(self.folderlist):
      f.write(fnames[i])
  @property
  def samples(self):
    return [f.samples for f in self.folderlist]
  def plot_getdist(self,ax=None,**kwargs):
    from getdist.plots import get_subplot_plotter
    gdfolders = [f.to_getdist() for f in self.folderlist]
    spp = get_subplot_plotter()
    spp.triangle_plot(gdfolders,filled=True,**kwargs)
  def to_getdist(self):
    return [f.to_getdist() for f in self.folderlist]
