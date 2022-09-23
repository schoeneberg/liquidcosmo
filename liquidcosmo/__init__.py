from .matplotlib_defaults import initialize_plots
initialize_plots()
from .folder import folder
from .foldercollection import foldercollection
def load(*args,**kwargs):
  if len(args)==1 and isinstance(args[0],str):
    return folder.load(args[0],kwargs)
  else:
    return foldercollection.load(*args,**kwargs)
def load_chain(path,**kwargs):
  return folder.load(path,**kwargs).get_chain()
def load_manually(paths,**kwargs):
  return folder.load_manually(path,**kwargs)
#variational NN similar to GP
