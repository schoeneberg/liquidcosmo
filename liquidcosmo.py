from src.folder import folder
def load(path):
  return folder.load(path)
def load_chain(path):
  return folder.load(path).get_chain()
#variational NN similar to GP
