from src.folder import folder
from src.matplotlib_defaults import matplotlib_defaults
matplotlib_defaults()
def load(path):
  return folder.load(path)
def load_chain(path):
  return folder.load(path).get_chain()
#variational NN similar to GP
