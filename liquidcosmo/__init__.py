from .matplotlib_defaults import initialize_plots
initialize_plots()
from .folder import folder
from .foldercollection import foldercollection

# Load a chain generally
def load(*args,**kwargs):
  if len(args)==1 and isinstance(args[0],str):
    return folder.load(args[0],**kwargs)
  else:
    return foldercollection.load(*args,**kwargs)

# Directly load chain only
def load_chain(path,**kwargs):
  return folder.load(path,**kwargs).get_chain()

# Load a file manually (providing all necessary keyword arguments)
def load_manually(paths,**kwargs):
  return folder.load_manually(path,**kwargs)

# Load a bestfit file from montepython as a single-item chain
def loadbestfit(path,**kwargs):
  return folder.loadbestfit(path,**kwargs)
#variational NN similar to GP

# TODO :: Update this method
def add_contours(getdist_plotter_instance, all_parnames, paramname, mean, sigmaupper = None, sigmalower = None, sigma = None, label=None, alpha1=0.3, alpha2=0.2, color="grey"):
  spp = getdist_plotter_instance
  if (sigmalower or sigmaupper) and sigma:
    raise Exception("Cannot pass 'sigma' and either of 'sigmalower' or 'sigmaupper'")
  elif sigma:
    sigmalower = sigma
    sigmaupper = sigma
  high = mean + sigmaupper
  low = mean-sigmalower
  newmean = 0.5*(high+low)
  newsig = 0.5*(high-low)

  for i in range(len(all_parnames)):
    index = all_parnames.index(paramname)
    if i < index:
        spp.add_y_bands(newmean, newsig, color=color, ax=spp.subplots[index,i], alpha1=alpha1, alpha2=alpha2, label = label,zorder=-2)
    if i > index:
        spp.add_x_bands(newmean, newsig, color=color, ax=spp.subplots[i,index], alpha1=alpha1, alpha2=alpha2, label = label,zorder=-2)

