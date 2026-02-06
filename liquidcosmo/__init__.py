from .settings import initialize_plots
initialize_plots()
from .folder import folder
from .foldercollection import foldercollection

from .settings import default_colors as colors
from .settings import initialize_plots as set_plot_settings

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

# Load from a cobaya sample object, or a pandas data object
def load_from(dataobject,**kwargs):
  return folder.load_from(dataobject,**kwargs)


#variational NN similar to GP

# TODO :: Update this method
def add_contours(getdist_plotter_instance, all_parnames, paramname, mean, sigmaupper = None, sigmalower = None, sigma = None, label=None, alpha1=0.3, alpha2=0.2, color="grey"):
  spp = getdist_plotter_instance
  if (sigmalower or sigmaupper) and sigma:
    raise ValueError("Cannot pass 'sigma' and either of 'sigmalower' or 'sigmaupper'")
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

def get_gaussian_chain(mean, *, std=None, cov=None, names=None, N=10000, bounds=None, **kwargs):
  import numpy as np
  mean = np.atleast_1d(mean)
  dim = mean.shape[0]
  if std is not None and cov is not None:
    raise ValueError("Cannot pass both 'std' and 'cov' at the same time!")
  elif cov is not None:
    cov = np.atleast_2d(cov)
    if cov.shape[0]!=cov.shape[1]:
      raise ValueError("The covariance matrix 'cov' needs to be square shaped!")
    if cov.shape[0]!=dim:
      raise ValueError("The covariance matrix 'cov' needs to have the same dimension as the mean vector 'mean'")
  elif std is not None:
    std = np.atleast_1d(std)
    if std.shape[0] !=dim:
      raise ValueError("The standard deviation vector 'std' needs to have the same dimension as the mean vector 'mean'")
    if std.size == 1:
      cov = np.eye(dim) * std[0]**2 # Assume same std for all items
    else:
      cov = np.diag(std**2)
  else:
    raise ValueError("Please pass either 'std' or 'cov'")
  samps = np.random.multivariate_normal(mean, cov, size=N)
  fo = load_from(samps, names=names,**kwargs)
  fo.chain._d['lnp'] = 0.5*np.sum((mean-samps) @ np.linalg.inv(cov) * (mean-samps), axis=1)
  if bounds is not None:
    bounds = np.atleast_2d(bounds)
    if bounds.shape[0]!=mean.shape[0]:
      raise ValueError("Please pass as many bounds as items in the mean")
    for name, bound in zip(fo.names[2:], bounds):
      fo.set_range(name, bound)
  return fo
