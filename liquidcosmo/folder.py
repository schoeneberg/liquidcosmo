import os
import numpy as np
import signal
import multiprocessing
from .chain import chain
import ast
from datetime import date
from collections import OrderedDict
from copy import deepcopy
from functools import partial

class _lq_code_type:
  montepython = 0
  cobaya = 1
  emcee = 2

class obj_iterator:
  def __init__(self,obj):
    self.current = 0
    self.N = obj.N
    self.obj = obj
  def __next__(self):
    if self.current<self.N:
      res = self.obj[self.current]
      self.current+=1
      return res
    else:
      raise StopIteration


def tex_convert(constr,withdollar,mean,texname=None,name=None,equalize=2.5):
  from .util import round_reasonable
  string = ""
  if constr[1]=="unconstrained":
    string= (texname+" " if texname is not None else "") +"unconstrained"
  elif constr[1]==">":
    string= (texname+" " if texname is not None else "") +"> "+round_reasonable(constr[0][0])
  elif constr[1]=="<":
    string= (texname+" " if texname is not None else "") +"< "+round_reasonable(constr[0][1])
  elif constr[1]=="+-":
    string= (texname+" = " if texname is not None else "") +round_reasonable(mean,errp=constr[0][1]-mean,errm=mean-constr[0][0],equalize=equalize)
  if withdollar:
    string = "$"+string+"$"
  if name is not None:
    string = name+":"+string
  return string

#from collections import OrderedDict
class folder:

  __limit_safety_factor = 0.7
  __equalize_errors = 2.5
  __precision_mode = False

  def __init__(self):
    self.verbose = 0
    self._foldername = None
    self._allchains = None
    self._chainprefix = None
    self._code = None
    self.lens = None
    self._texnames = None
    self._arr = None
    self._narr = None
    self._log =  None
    self._confinfo = {}
    self.path = None
    self.tag = None

  # copy this object
  def copy(self):
    a = folder()
    a.verbose = self.verbose
    a._foldername = self._foldername
    a._allchains = self._allchains
    a._chainprefix = self._chainprefix
    a._code = self._code
    a.lens = self.lens
    a._texnames = self._texnames
    a._arr = self._arr
    a._narr = self._narr
    a._log = self._log
    a._confinfo = self._confinfo
    a.path = self.path
    a.tag = self.tag
    return a

  # deepcopy this object
  def deepcopy(self):
    a = folder()
    a.verbose = self.verbose
    a._foldername = deepcopy(self._foldername)
    a._allchains = deepcopy(self._allchains)
    a._chainprefix = deepcopy(self._chainprefix)
    a._code = deepcopy(self._code)
    a.lens = deepcopy(self.lens)
    a._texnames = deepcopy(self._texnames)
    a._arr = deepcopy(self._arr)
    a._narr = deepcopy(self._narr)
    a._log = deepcopy(self._log)
    a._confinfo = deepcopy(self._confinfo)
    a.path = deepcopy(self.path)
    a.tag = deepcopy(self.tag)
    return a

  # Construct a folder object (containing multiple physical chains) from a path
  @classmethod
  def load(obj, path, kind="all", burnin_threshold = 3, verbose = 0,timeout=60, tag=None, keep_non_markovian=True):
    a = obj()
    a.verbose = verbose
    a._foldername, a._allchains, a._chainprefix, a._code = obj._resolve_chainname(path, kind=kind)
    a.lens = None
    a._texnames = None
    a._arr = None
    a._narr = None
    a._log = None
    a._confinfo = {}
    a.path = os.path.abspath(a._foldername)
    if tag==None:
      if a._code==_lq_code_type.montepython:
        a.tag = os.path.basename(os.path.dirname(a.path) if not os.path.isdir(a.path) else a.path)
      elif a._code==_lq_code_type.cobaya:
        a.tag = os.path.basename(a._chainprefix)
      elif a._code==_lq_code_type.emcee:
        a.tag = os.path.basename(a._chainprefix)
    else:
      a.tag = tag
    a.get_chain(burnin_threshold=burnin_threshold,timeout=timeout,keep_non_markovian=keep_non_markovian)
    return a

  @classmethod
  def loadbestfit(obj, path, verbose=0, tag=None):
    a = obj()
    a.verbose = verbose
    if not os.path.isfile(path):
      raise Exception("Need to point 'loadbestfit' function to a bestfit file. Could not open {}".format(path))
    a._foldername = os.path.dirname(path)
    a._allchains = []
    a._chainprefix = ""
    a._code = _lq_code_type.montepython
    a.lens = [1]
    a._arr = None
    a._narr = None
    a._log = None
    a._confinfo = {}
    a.path = os.path.abspath(a._foldername)
    if tag==None:
      if a._code==_lq_code_type.montepython:
        a.tag = os.path.basename(os.path.dirname(a.path) if not os.path.isdir(a.path) else a.path)
      elif a._code==_lq_code_type.cobaya:
        a.tag = os.path.basename(a._chainprefix)
      elif a._code==_lq_code_type.emcee:
        a.tag = os.path.basename(a._chainprefix)
    else:
      a.tag = tag
    arrdict = OrderedDict({'N':np.array([1],dtype=int),'lnp':np.array([0],dtype=float)})
    with open(path,"r") as inf:
      names = [x.strip() for x in inf.readline()[1:].split(",")]
      vals = [float(x) for x in inf.readline().split() if x]
      for name, val in zip(names,vals):
        arrdict[name] = np.array([val],dtype=float)
    a._narr = chain(arrdict)
    a._texnames = names
    return a

  @classmethod
  def load_manually(obj, filenames,verbose=0, names=None,**kwargs):
    if isinstance(filenames,str):
      filenames = [filenames]
    a = obj()
    a.verbose = verbose
    a._foldername = "MANUALLY_LOADED"
    a._allchains = filenames
    a._chainprefix = ""
    a.lens = None
    a._texnames = None
    a._arr = None
    a._narr = None
    a._log = None
    a.path = ""
    arr = a._get_array(**kwargs)
    arrprime = arr
    if names:
      if not isinstance(names[0],str):
        raise Exception("The argument 'names' has to be a list of strings")
      if not "N" in names:
        raise Exception("The argument 'names' has to be a list of strings, containing at least 'N' and 'loglkl' or 'chi2' (missing: 'N')")
      if (not "loglkl" in names) and (not "chi2" in names):
        raise Exception("The argument 'names' has to be a list of strings, containing at least 'N' and 'loglkl' or 'chi2' (missing: 'loglkl' or 'chi2')")
      if "loglkl" in names:
        idx_loglkl = np.argmax(names=="loglkl")
        fac_loglkl = 1.
      else:
        idx_loglkl = np.argmax(names=="chi2")
        fac_loglkl = 0.5
      idx_N = np.argmax(names=="N")
    else:
      idx_loglkl = 1
      idx_N = 0
      fac_loglkl = 1.
      names = []

    arrdict = OrderedDict({'N':arr[idx_N],'lnp':arr[idx_loglkl]*fac_loglkl})
    for i in range(len(arr)):
      #Fill in known parameters
      if i<len(names):
        if i==idx_loglkl or i==idx_N:
          continue
        arrdict[names[i]]=arr[i]
      #Fill in unknown parameters
      else:
        arrdict["UNKNOWN.{}".format(i-len(names))]=arr[i]
    a._narr = chain(arrdict)
    a._texnames = {name:name for name in arrdict.keys()} #we don't try any re-naming
    a._log = {}
    return a

  # Construct a folder object (containing multiple physical chains) from a data object
  @classmethod
  def load_from(obj, dataobj, verbose=0, burnin_threshold = 3, tag=None):
    a = obj()
    a.verbose = verbose
    a._foldername, a._allchains, a._chainprefix, a._code = "from data object",["from data object"],"from data object",_lq_code_type.montepython
    a.path = "from data object"
    a.tag = tag if tag else "from data object"
    a.convert_chain(dataobj,burnin_threshold=burnin_threshold)
    return a

  # -- Load a given chain (for a given full filename)
  def __ldchain__(filename,verbose=False,precision_mode=False,checkbroken=False,keep_non_markovian=True):
    # This should all be effectively equivalent to (but vastly faster and more memory efficient than)
    #arr = np.genfromtxt(filename).T
    if verbose:
      print("liquidcosmo :: Loading chain {}".format(filename))

    # Fast reading of the file to get line count
    def fastread(f):
      f.seek(0)
      # Helper function to put a large file into blocks of decent size
      def blocks(f, size=65536): #2^16
        while True:
          b = f.read(size)
          if not b: break
          yield b
      # use the sum over a generator in order to avoid loading all at once
      return sum(bl.count('\n') for bl in blocks(f))

    # Reading of first line of elements to get number of elems per line
    def get_num_items(f):
      f.seek(0)
      for line in f:
        if line.startswith("#"):
          continue
        else:
          return len(line.split())

    # Putting the actual contents of the file into an array called 'arr' (HAS TO EXIST!!)
    def load_file(f,storage_type):
      f.seek(0)
      i = 0
      nm = 0
      for line in f:
        if line.startswith('#'):
          if 'update proposal' in line:
            nm = i
          continue
        arr[i] = [storage_type(f) for f in line.split()]
        i+=1
      return i, nm

    # A check to see if the last line of a file is broken for some reason
    def check_last_line(f,numitems):
      check_last_line = f.seek(0, 2)
      current = f.tell()
      #should be MORE than enough to catch last line (40 characters per item)
      leng = numitems*40
      if current < leng:
        leng=current
      f.seek(current-leng)
      read = f.read(leng)
      if read.count('\n')<1:
        return False
      else:
        return len(read.split("\n")[-1].split())==0 and len(read.split("\n")[-2].split())==numitems

    # The actual code of the function
    with open(filename,"r",encoding="utf-8",errors='ignore') as f:
      linenum = fastread(f)
      numitems = get_num_items(f)
      if linenum<=0 or numitems==0 or numitems==None:
        return []
      if checkbroken and not check_last_line(f,numitems):
        raise Exception("File broken :: {} (last line is not proper)".format(filename))
      # The chain files are usually written at a relatively low precision, allowing us to use low precision representations. However, we still give the user the option to use full precision
      storage_type = (np.float64 if precision_mode == True else np.float32)
      arr = np.empty((linenum,numitems),dtype=storage_type)
      real_len, non_markovian = load_file(f,storage_type) # MODIFIES arr (!!)
      if keep_non_markovian:
        arr = arr[:real_len]
      else:
        arr = arr[non_markovian:real_len]

    # Just some verbosity
    if verbose:
      print("liquidcosmo :: Success in loading chain {}".format(filename))

    return arr.T



  # -- Resolve the chain names of chains --
  @classmethod
  def _resolve_chainname(obj,path,kind="all"):
    if os.path.isfile(path):
      folder = os.path.dirname(path)
      allchains = [path]
      code = _lq_code_type.montepython
    elif os.path.isdir(os.path.join("chains",path)):
      folder = os.path.join("chains",path)
      allchains, code = obj._get_chainnames(folder,kind=kind)
    elif os.path.isdir(path):
      folder = path
      allchains, code = obj._get_chainnames(folder,kind=kind)
    elif os.path.exists(os.path.join("chains",path)):
      folder = os.path.dirname(os.path.join("chains",path))
      allchains, code = obj._get_chainnames(folder,kind=kind)
    elif os.path.exists(os.path.join("chains",path+"__1.txt")):
      folder = os.path.dirname(os.path.join("chains",path+"__1.txt"))
      allchains, code = obj._get_chainnames(folder,kind=kind)
      allchains = [ch for ch in allchains if path in ch]
    else:
      folder = os.path.dirname(path)
      if folder=="" or not os.path.isdir(folder):
        raise Exception("Could not find a chain from (check the folder name): "+path)
      allchains, code = obj._get_chainnames(folder,kind="match:"+os.path.basename(path))
    prefix = obj._resolve_prefix(np.sort(allchains)[0])
    return folder,allchains,prefix,code

  @classmethod
  def _resolve_prefix(obj, path):
    if "__" in path:
      return path.split("__")[0]+"_"
    elif "_1.txt" in path:
      return path.split("_1.txt")[0]
    elif ".txt" in path:
      return ".".join(path.split(".")[:-2])
    elif ".hdf5" in path:
      return path.split(".hdf5")[0]
    else:
      raise Exception("Unrecognized path for prefix recognition (please report to developer): "+str(path))
  # -- Get the names of the given chains
  @classmethod
  def _get_chainnames(obj,directory,kind="newest"):

    # Check all files in directory
    allfilenames  = os.listdir(directory)

    # Detect the code, for now only montepython and cobaya are supported
    if np.any(["__" in chain for chain in allfilenames]):
      code = _lq_code_type.montepython
      chain_targets = [chain for chain in allfilenames if "__" in chain]
      delim = "__"
    elif np.any(["_1.txt" in chain for chain in allfilenames]):
      code = _lq_code_type.montepython
      chain_targets = [chain for chain in allfilenames if np.any(["_{}.txt".format(i) in chain for i in range(100)])]
      delim = "_"
    elif np.any([".hdf5" in chain for chain in allfilenames]):
      code = _lq_code_type.emcee
      chain_targets = [chain for chain in allfilenames if ".hdf5" in chain]
      delim = ".hdf5"
    else:
      code = _lq_code_type.cobaya
      chain_targets = [chain for chain in allfilenames if ".txt" in chain]

    if kind=="newest" and code== _lq_code_type.montepython:
      chain_targets = np.sort(chain_targets)[::-1]
      newest_chain_dates = chain_targets[0].split("__")[0]
      allchains= [os.path.join(directory,chain) for chain in chain_targets if newest_chain_dates in chain]

    elif kind=="all":
      allchains= [os.path.join(directory,chain) for chain in chain_targets]
    elif kind.startswith("match:"):
      allchains= [os.path.join(directory,chain) for chain in chain_targets if kind.replace("match:","") in chain]

    else:
      raise Exception("Kind can either be 'all' or 'newest'")

    return allchains, code

  # -- Load all points from folder --
  def _get_array(self,excludesmall=True,burnin_threshold=3,timeout=60,keep_non_markovian=True):
    if self._arr is not None:
      return self._arr

    chainnames = self._allchains

    if excludesmall:
      chainnames = [chain for chain in chainnames if not "_1__" in chain]

    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(8)
    signal.signal(signal.SIGINT, sigint_handler)
    loading_function = partial(folder.__ldchain__,precision_mode=self.__precision_mode,keep_non_markovian=keep_non_markovian)
    try:
      filearr = pool.map_async(loading_function, chainnames)
      filearr = filearr.get(timeout) # by default 60 seconds (timeout)
    except KeyboardInterrupt:
      print("You stopped the loading of the chains by KeyboardInterrupt")
      pool.terminate()
      raise
    else:
      pool.close()
    pool.join()

    filearr = [fa for fa in filearr if not ((isinstance(fa,list) and fa==[]) or fa.ndim<=1)]

    #filearr = [fa for fa in filearr if fa!=[] and fa.ndim>1]
    if(len(filearr))==0:
      raise Exception("There is probably a problem with the chain folder '{}' that is attempted to be analyzed. Please make sure that the folder is non-empty and the chains are not empty files.".format(self._foldername))

    if burnin_threshold >= 0:
      import scipy.stats
      total_removed = 0
      total_len = 0
      #thres = scipy.stats.chi2.isf(scipy.special.erfc(burnin_threshold/np.sqrt(2)),len(arrs)-2)
      for j in range(len(filearr)):
        idx = np.argmax(filearr[j][1]<np.min(filearr[j][1])+burnin_threshold)
        filearr[j] = filearr[j][:,idx:]
        total_removed += idx
        total_len += len(filearr[j][0])

      if self.verbose > 0:
        print("liquidcosmo :: Removed burnin [%i/%i] for chain in folder '%s'."%(total_removed,total_len,self._foldername))

    self.lens = np.array([len(fa[0]) for fa in filearr])


    # The code below should be completely equivalent to
    #arrs = np.hstack(filearr)
    # but that causes memory problems, since it isn't smart enough to lazily allocate
    # (The problem cannot be avoided in its entirety of course)
    arrs = np.empty((len(filearr[0]),np.sum(self.lens)))

    idx = 0
    # This is technically not the most efficient loop, but it avoids memory issues
    for j in range(len(filearr)):
      leng = self.lens[j]
      for iparam in range(len(filearr[j])):
        arrs[iparam][idx:idx+leng] = filearr[j][iparam][:]
      idx+=leng
      # Important: Free memory that we don't need anymore!
      filearr[j] = None

    self._arr = arrs
    return arrs

  # -- Get dictionary of all parameters in .paramnames file --
  # Returns a dictionary with parameter name and index
  def _load_names(self):
    retdict = {}
    texdict = {'N':'Multiplicity','lnp':'-\\ln(\\mathcal{L})'}
    index = 2
    if self._code == _lq_code_type.montepython:
      with open(self._chainprefix+".paramnames") as parnames:
        line = parnames.readline()
        while line:
          texname = " ".join(line.split()[1:]).strip()
          paramname = line.split()[0].strip()
          if self.verbose>2:
            print("liquidcosmo :: Param {} was found at index {}".format(paramname,index))
          retdict[paramname]=index
          texdict[paramname] = texname
          index+=1
          line = parnames.readline()
    elif self._code == _lq_code_type.cobaya:
      with open(self._chainprefix+".1.txt") as parnames:
        line = parnames.readline()
        if line[0]!="#":
          raise Exception("Malformed cobaya chain file? No header detected!")
        names = [x for x in line[1:-1].split(" ") if x]
        retdict = {name:i for i,name in enumerate(names)}
        texdict = {name:name for name in names}
        if self.logfile != {}:
          texdict = {name:self.logfile['parinfo'].get(name,{}).get('latex',name) for name in names}
    else:
      raise Exception("Missing code type (report to developer)")
    return retdict,texdict

  def get_chain_emcee(self):
    if not (self._narr is None):
      return self._narr
    all_samples = []
    all_log_prob_samples = []
    for fname in self._allchains:
      import emcee
      reader = emcee.backends.HDFBackend(fname)
      tau = reader.get_autocorr_time(tol=0)

      if not any(np.isnan(tau)):
        burnin = int(2 * np.max(tau))
        thin = int(0.5 * np.min(tau))
      else:
        burnin =  int(0.2*reader.iteration)

      thin = 1
      samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
      log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)

      mask = log_prob_samples > np.max(log_prob_samples)-100
      all_samples.append(samples[mask,:])
      all_log_prob_samples.append(log_prob_samples[mask])
    samples = np.vstack(all_samples)
    log_prob_samples = np.concatenate(all_log_prob_samples)
    arrdict = OrderedDict({'N':np.ones(len(log_prob_samples)),'lnp':log_prob_samples})
    for ip in range(len(samples.T)):
      arrdict['x{}'.format(ip)] = samples[:,ip]
    self._narr = chain(arrdict)
    self._texnames = {name:name for name in arrdict}

  # Create the chain object (equivalent to a named dictionary or arrays)
  def get_chain(self,excludesmall=True,burnin_threshold=5,timeout=60,keep_non_markovian=True):
    if not (self._narr is None):
      return self._narr
    if self._code == _lq_code_type.emcee:
      return self.get_chain_emcee()
    arr = self._get_array(excludesmall=excludesmall,burnin_threshold=burnin_threshold,timeout=timeout,keep_non_markovian=keep_non_markovian)
    arrdict = OrderedDict({'N':arr[0],'lnp':arr[1]})
    #arrdict = {'N':arr[0],'lnp':arr[1]}
    parnames, texnames = self._load_names()
    for index in range(2,len(arr)):
      found = False
      for key,value in parnames.items():
        if value==index:
          if not found:
            arrdict[key]=arr[index]
            found=True
          else:
            raise Exception("Key {} was found twice (already in {})".format(key,arrdict.keys()))
      if not found:
        arrdict["UNKNOWN.{}".format(index)]=arr[index]
    self._narr = chain(arrdict)
    self._texnames = texnames
    return self._narr

  # Create the chain object (equivalent to a named dictionary or arrays)
  def convert_chain(self,dataobj,burnin_threshold=5):
    if not (self._narr == None):
      return self._narr

    if hasattr(dataobj,"_data") and not hasattr(dataobj, "_get_numeric_data"):
      data = getattr(dataobj,"_data")
    else:
      data = dataobj

    names = list(data.columns)

    if (self._arr == None):

      self.lens = [len(data.index)]

      self._arr = [data[name].to_numpy() for name in names]

    parnames = {name:i for i,name in enumerate(names)}

    arr = self._arr

    if (names[0]!='N' and names[0]!='weight'):
      raise ValueError("Needs either 'N' or 'weight' in the first position of the loaded object, but only has '{}'".format(names))
    if (names[1]!='lnp' and names[1]!='logp' and names[1]!='loglike' and names[1]!='minuslogpost'):
      raise ValueError("Needs any of ['lnp','logp','loglike','minuslogpost'] in the second position of the loaded object, but only has '{}'".format(names))

    arrdict = OrderedDict({'N':arr[0],'lnp':arr[1]})
    for index in range(2,len(arr)):
      found = False
      for key,value in parnames.items():
        if value==index:
          if not found:
            arrdict[key]=arr[index]
            found=True
          else:
            raise KeyError("Key {} was found twice (already in {})".format(key,arrdict.keys()))
      if not found:
        arrdict["UNKNOWN.{}".format(index)]=arr[index]
    self._narr = chain(arrdict)
    self._texnames = {name:name for name in names}

    return self._narr

  # -- Whatever you do to a folder, typically you want to do to the underlying chain
  def __getitem__(self,q):
    if isinstance(q,(int,np.integer)):
      return self.chain.get_dict(q)
    elif not isinstance(q,str):
      res = self.copy()
      res._narr = self.chain[q]
      return res
    else:
      if q not in self.names:
        raise KeyError("Key '{}' not found within the chain".format(q))
      return self.chain[q]
  def __setitem__(self,q,v):
    empty = {'log':0,'initial':1,'bound':[None,None],'initialsigma':1,'type':'derived'}
    if isinstance(q,str):
      if not q in self._texnames:
        self._texnames[q] = q
        if self.logfile != {}:
          if q in self.logfile['parinfo']:
            if not self.logfile['parinfo'][q]['initialsigma']==0:
              raise Exception("Did not expect parameter '{}' in logfile. This is a bug, please report to the developer".format(q))
          else:
            self.logfile['parinfo'][q] = empty
      elif (self.logfile!={}) and (q not in self.logfile['parinfo']):
        self.logfile['parinfo'][q] = empty
    self.chain[q] = v

  def __contains__(self,m):
    return q in self.names

  @property
  def bestfit(self):
    res = self.copy()
    res._narr = self.chain.bestfit
    return res
  def __str__(self):
    return "Folder"+self.chain._str_part()
  def __repr__(self):
      return self.__str__()

  # Basically same as setting a value directly, but also allows passing a function
  def derive(self, name, func, texname = None):
    if self.verbose > 1:
      print("liquidcosmo :: Deriving new parameter "+name)
    self._texnames[name] = (name if texname is None else texname)
    if isinstance(func,(list,tuple,np.ndarray)):
      self[name] = func
    self[name] = func(self)
  @property
  def cosmopars(self):
    if self.logfile == {}:
      return None
    else:
      return [x for (x,v) in self.logfile['parinfo'].items() if v['type']=='cosmo']
  @property
  def cosmoargs(self):
    if self.logfile == {}:
      return None
    else:
      return list(self.logfile['arginfo'].keys())
  def __iter__(self):
    return obj_iterator(self)
  @property
  def N(self):
    return self.chain.N
  @property
  def logfile(self):
    if self._log is None:
      self._log = self._read_log()
      return self._log
    else:
      return self._log
  @property
  def names(self):
    return self.chain.names
  @property
  def chain(self):
    return self.get_chain()
  @property
  def samples(self):
    return np.array([self.chain._d[arg] for arg in self.names[2:]])
  @property
  def d(self):
    return len(self.names)-2

  def _read_log(self):
    if self._code == _lq_code_type.montepython:
      return self._read_log_montepython()
    elif self._code == _lq_code_type.cobaya:
      return self._read_log_cobaya()
    elif self._code == _lq_code_type.emcee:
      return {}
    else:
      raise Exception("Unexpected code type")

  # -- Get all content of a log file capturing important information about the sampling process
  # -- Currently, only montepython log.param files are supported
  def _read_log_montepython(self):
    self._log = {}
    try:
      loginfo ={'path':{}}
      parinfo = {}
      arginfo = {}
      lklopts = {}
      with open(os.path.join(self.path,'log.param')) as logfile:
        line = logfile.readline()
        while line:
          line = line.strip()
          if("#-----CLASS" in line):
            before,after = line.split("(")
            after = after.split(")")[0].strip()
            loginfo["version"] = before.split("CLASS")[1].strip()
            one,two=after.split(",")
            one=one.split(":")
            two=two.split(":")
            if "branch" in one[0]:
              loginfo["branch"]=one[1].strip()
            if "hash" in two[0]:
              loginfo["hash"]=two[1].strip()
          elif line=="" or line[0]=='#':
            line = logfile.readline()
            continue
          elif "data.experiments" in line:
            loginfo["experiments"]=ast.literal_eval(line.split("=")[1])
          elif "data.over_sampling" in line:
            loginfo["oversampling"]=ast.literal_eval(line.split("=")[1])
          elif "data.parameters" in line:
            before,after = line.split("=")
            parname = ast.literal_eval(before.split("[")[1].split("]")[0])
            parinfo[parname] = {}
            afterlist = after.split("[")[1].split("]")[0].split(",")
            if len(afterlist)==6:
              for i in range(len(afterlist)):
                afterlist[i] = afterlist[i].strip()
              parinfo[parname]["initial"]=float(afterlist[0])*float(afterlist[4])
              parinfo[parname]["bound"]=[ast.literal_eval(afterlist[1]),ast.literal_eval(afterlist[2])]
              #if parinfo[parname]["bound"][0] is not None:
              #  parinfo[parname]["bound"][0]*=float(afterlist[4])
              #if parinfo[parname]["bound"][1] is not None:
              #  parinfo[parname]["bound"][1]*=float(afterlist[4])
              parinfo[parname]["initialsigma"]=float(afterlist[3])*float(afterlist[4])
              parinfo[parname]["type"]=ast.literal_eval(afterlist[5])
            else:
              raise Exception("Unknown argument list in log.param line:: "+repr(line))
            if parname[:3]=="log":
              parinfo[parname]['log']=10
            elif parname[:2]=="ln":
              parinfo[parname]['log']=np.e
            else:
              parinfo[parname]['log']=0
          elif "data.cosmo_arguments" in line and "data.path['cosmo']" in line:
            pass
          elif "data.cosmo_arguments" in line and not "update" in line:
            before,after = line.split("=")
            parname = ast.literal_eval(before.split("[")[1].split("]")[0])
            arginfo[parname] = {}
            arginfo[parname]=ast.literal_eval(after.strip())
          elif "data.cosmo_arguments.update" in line:
            loginfo["command"] = ast.literal_eval(line.split("(",1)[1].rsplit(")",1)[0])
          elif "data.emulator" in line:
            if not "emulator" in loginfo:
              loginfo["emulator"] = {}
            before, after = line.split("=")
            emuinfo_par = ast.literal_eval(before.split("[")[1].split("]")[0])
            loginfo["emulator"][emuinfo_par] = after.strip()
          elif "data.path" in line:
            before,after = line.split("=")
            pathname = before.split("[")[1].split("]")[0]
            loginfo['path'][ast.literal_eval(pathname)] = ast.literal_eval(after.strip())
          elif "=" in line:
            before,after = line.split("=")
            lkl,optname = before.strip().split(".")
            try:
              lklopts[lkl]
            except:
              lklopts[lkl] = {}
            if after.strip()=="Ellipsis":
              after = "..."
            lklopts[lkl][optname] = ast.literal_eval(after.strip())
          else:
            raise Exception("Unrecognized line in log.param :\n",repr(line))
          line = logfile.readline()
      if self.verbose>2:
        print("loginfo = ",loginfo)
        print("parinfo = ",parinfo)
        print("arginfo = ",arginfo)
        print("lklopts = ",lklopts)
      self._log = {'loginfo':loginfo,'parinfo':parinfo,'arginfo':arginfo,'lklopts':lklopts}
    except Exception as e:
      print("ERROR reading logfile : ",e)
      self._log = {}
    return self._log

  def _read_log_cobaya(self):
    self._log = {}
    loginfo ={'path':{}}
    parinfo = {}
    arginfo = {}
    lklopts = {}
    with open(self._chainprefix+'.updated.yaml') as logfile:
      try:
        import yaml
      except ImportError as ie:
        raise Exception("Currently cannot run without pyyaml, sorry about that! Please run 'pip install pyyaml'") from ie
      logdict = yaml.safe_load(logfile)
      lklopts = logdict.pop('likelihood')
      parinfo = logdict.pop('params')
      for pname in parinfo:
        par = parinfo[pname]
        par['bound'] = [None,None]
        if 'prior' in par and 'dist' in par['prior'] and par['prior']['dist'] == 'uniform':
          par['min'] = par['prior'].pop('loc',par['prior'].pop('min',0))
          par['max'] = par['min']+par['prior'].pop('scale',par['prior'].pop('max',1))
        if 'min' in par:
          par['bound'][0]  = par['min']
        if 'max' in par:
          par['bound'][1]  = par['max']
      arginfo = logdict.pop('theory')
      loginfo = logdict
      self._log = {'loginfo':loginfo,'parinfo':parinfo,'arginfo':arginfo,'lklopts':lklopts}
    return self._log

  def __get_individual_counts(self):
    if self.N == np.sum(self.lens):
      # If we still have the original number of points, the answer is trivial
      counts = self.lens
    else:
      # If the chain has been 'tampered' with, we need to re-distribute the remaining point over the given number of files
      if(self.N < len(self.lens)):
        # Less chain points than files, summarize all in one file
        counts = [self.N]
      else:
        counts = np.empty_like(self.lens)
        c_per_file = self.N//len(self.lens)
        remainder = self.N-c_per_file*len(self.lens)
        for i in range(len(self.lens)):
          counts[i] = c_per_file + (1 if remainder>0 else 0)
          if remainder > 0:
            remainder -= 1
    return counts

  def to_individual(self):
    counts = self.__get_individual_counts()

    from .foldercollection import foldercollection
    coll = foldercollection()
    index = 0
    for i,c in enumerate(counts):
      fo = self.deepcopy() #copy all info (like log.param), rest will be changed below anyway
      chain = self.chain[index:index+c]
      fo._narr = chain
      fo._arr = None #Technically not necessary, but let's do it for safety
      fo.lens = np.array([c])
      fo._allchains = [self._allchains[i]]
      coll.folderlist.append(fo)
      index+=c
    return coll

  def cut(self,first=0.0,last=1.0,thin=1):
    if len(self.lens)==1:
      fo = self.deepcopy()
      fo = fo[int(fo.N*first):int(fo.N*last)+1:thin]
      fo.lens[0] = fo.N
      if(fo.N<1):
        raise ValueError("Could not cut array, since less than 1 element survived")
      return fo
    else:
      fc = self.deepcopy().to_individual()
      for i in range(len(fc.folderlist)):
        fc.folderlist[i] = fc.folderlist[i].cut(first=first,last=last,thin=thin)
      return fc.merge()

  def write(self, fname, codetype="montepython"):
    if not (codetype=="montepython" or codetype=="Montepython" or codetype=="MontePython" or codetype=="MONTEPYTHON"):
      raise Exception("Can currently only write in MontePython style")
    if self._code != _lq_code_type.montepython:
      raise Exception("Sadly this conversion from code type '{}' to 'montepython' is not yet supported".format(self._code))
    if os.path.exists(fname):
      raise Exception("File path already exists : ",fname)
    else:
      # How many points to put into the individual chains?
      counts = self.__get_individual_counts()
      os.mkdir(fname)
      index = 0
      for i, c in enumerate(counts):
        with open(os.path.join(fname,str(date.today())+"_"+str(c)+"__"+str(i+1)+".txt"),"w") as ofile:
          arr = self.chain[index:index+c+1]._d
          index+=c
          ofile.write("# "+"\t".join(arr.keys())+"\n")
          for j in range(len(arr['N'])):
            ofile.write("%.4g"%arr['N'][j])
            ofile.write(" %.6g\t"%arr['lnp'][j])
            for k in list(arr.keys())[2:]:
              ofile.write("%.6e\t"%arr[k][j])
            ofile.write("\n")
      with open(os.path.join(fname,str(date.today())+"_"+str(counts[0])+"_.paramnames"),"w") as ofile:
        for k in self.names[2:]:
          ofile.write("{} {}\n".format(k,self._texnames[k]))
      if self.logfile != {}:
        loginfo,parinfo,arginfo,lklopts = self.logfile['loginfo'],self.logfile['parinfo'],self.logfile['arginfo'],self.logfile['lklopts']
        with open(os.path.join(fname,'log.param'),"w") as logfile:
          logfile.write("#-----CLASS {} (branch: {}, hash: {})-----\n\n".format(loginfo['version'],loginfo['branch'],loginfo['hash']))
          logfile.write("data.experiments=[{}]\n\n".format(",".join("'"+str(x)+"'" for x in loginfo['experiments'])))
          logfile.write("data.over_sampling=[{}]\n\n".format(",".join(str(x) for x in loginfo['oversampling'])))
          for par in parinfo.keys():
            logfile.write("data.parameters['{}'] = [{}, {}, {}, {}, {}, '{}']\n".format(par,parinfo[par]['initial'],parinfo[par]['bound'][0],parinfo[par]['bound'][1],parinfo[par]['initialsigma'],1,parinfo[par]['type']))
          logfile.write("\n")
          for par in arginfo.keys():
            logfile.write("data.cosmo_arguments['{}'] = '{}'\n".format(par,arginfo[par]))
          logfile.write("\n")
          if 'command' in loginfo.keys():
            logfile.write("data.cosmo_arguments.update({})\n".format(str(loginfo['command'])))
          logfile.write("\n")
          for par in loginfo['path'].keys():
            logfile.write("data.path['{}'] = '{}'\n".format(par,loginfo['path'][par]))
          logfile.write("\n")
          for par in lklopts.keys():
            for k in lklopts[par]:
              lklop = lklopts[par][k]
              logfile.write("{}.{} = {}\n".format(par,k,str(lklop) if not isinstance(lklop,str) else "'{}'".format(lklop)))

  def set_range(self,parname,lower=None,upper=None,destructive=False):
    if isinstance(lower,list) or type(lower) is np.ndarray:
      if upper:
        raise ValueError("Cannot have 'upper!=None' but also 'lower' be a list")
      elif len(lower)!=2:
        raise ValueError("If you provide a list, it must be exactly 2 elements long")
      else:
        upper = lower[1]
        lower = lower[0]
    if lower is not None and destructive:
      self._narr = self.chain[self[parname]>=lower]
    if upper is not None and destructive:
      self._narr = self.chain[self[parname]<=upper]
    if self.logfile != {} and parname in self.logfile['parinfo']:
      self.logfile['parinfo'][parname]['bound'] = [lower,upper]
  def get_range(self,parname):
    if self.logfile != {} and parname in self.logfile['parinfo']:
      return self.logfile['parinfo'][parname]['bound']
    else:
      return [None,None]
  def range_mask(self,parname):
    lower,upper = self.get_range(parname)
    mask = np.ones(self.N,dtype=bool)
    if lower is not None:
      mask = np.logical_and(mask,self[parname]>=lower)
    if upper is not None:
      mask = np.logical_and(mask,self[parname]<=upper)
    return mask
  def cut_to_range(self, parname):
    mask = self.range_mask(parname)
    return self[mask]
  def get_bounds(self):
    bounds = {}
    for par in self.names[2:]:
      bounds[par] = self.get_range(par)
    return bounds

  def set_texname(self,parname,texname):
    if parname in self._texnames:
      self._texnames[parname] = texname
    else:
      raise Exception("Parameter '{}' not found in the list of parameters.".format(parname))

  def get_masked(self, mask):
    return self[mask].deepcopy()
  def subrange(self, parname=None, value=None):
    if isinstance(parname,str):
      if isinstance(value,(list,np.ndarray,tuple)):
        if len(value)!=2:
          raise Exception("Could not understand subrange arugment '{}' for parameter '{}'".format(value,parname))
        else:
          if isinstance(value[0],(float,np.floating)) and isinstance(value[1],(float,np.floating)):
            mask = np.logical_and(self[parname]>value[0],self[parname]<value[1])
            return self[mask]
          elif isinstance(value[0],(float,np.floating)) and value[1] is None:
            return self[self[parname]>value[0]]
          elif isinstance(value[1],(float,np.floating)) and value[0] is None:
            return self[self[parname]<value[1]]
          else:
            raise Exception("Could not understand subrange arugment '{}' for parameter '{}'".format(value,parname))
    elif isinstance(parname,(list,np.ndarray,tuple)):
      obj = self
      for i,p in enumerate(parname):
        obj = obj.subrange(p, value[i])
      return obj
    elif isinstance(parname,dict):
      obj = self
      for p,v in parname.items():
        obj = obj.subrange(p, v)
      return obj

  def mean(self,parnames=None,asdict=False):
    if asdict:
      mean = self.mean(parnames=parnames,asdict=False)
      return {p:mean[ip] for ip,p in enumerate(parnames or self.names[2:])}
    if isinstance(parnames,str):
      return self._parmean(parnames)
    elif isinstance(parnames,(list,tuple)):
      return np.array([self._parmean(q) for q in parnames])
    elif parnames is None:
      return np.array([self._parmean(q) for q in self.names[2:]])
    else:
      raise Exception("Input to mean '{}' not recognized.".format(parnames))

  def cov(self,parnames=None):
    if isinstance(parnames,str):
      return self._parcov([parnames])
    elif isinstance(parnames,(list,tuple)):
      return self._parcov(parnames)
    elif parnames is None:
      return self._parcov(self.names[2:])
    else:
      raise Exception("Input to mean '{}' not recognized.".format(parnames))
  def std(self,parnames=None,asdict=False):
    if asdict:
      std = self.std(parnames=parnames,asdict=False)
      return {p:std[ip] for ip,p in enumerate(parnames or self.names[2:])}
    if isinstance(parnames,str):
      return self._parstd(parnames)
    elif isinstance(parnames,(list,tuple)):
      return np.array([self._parstd(q) for q in parnames])
    elif parnames is None:
      return np.array([self._parstd(q) for q in self.names[2:]])
    else:
      raise Exception("Input to std '{}' not recognized.".format(parnames))

  def _parmean(self,parname):
    return np.average(self[parname],weights=self['N'])
  def _parcov(self,parnames):
    if self.N==1:
      raise Exception("Cannot compute covariance of a chain with only a single point!")
    return np.cov([self[parname] for parname in parnames],fweights=self['N'])
  def _parstd(self,parname):
    if self.N==1:
      raise Exception("Cannot compute covariance of a chain with only a single point!")
    mean = self._parmean(parname)
    return np.sqrt(np.average(np.array(self[parname]-mean)**2,weights=self['N']))

  def _upper_limit(self,parname, alpha):
    # By default, this computes a given quantile
    mask = self.range_mask(parname)
    Nmask = np.count_nonzero(mask)
    if parname in self._confinfo and Nmask == self._confinfo[parname]['Nmask']:
      samples = self[parname][mask]
      idxs = self._confinfo[parname]['idxs']
      wfracs = self._confinfo[parname]['wfracs']
    else:
      samples = self[parname][mask]
      idxs = np.argsort(samples)
      sorted_weights = self['N'][mask][idxs]
      wfracs = np.cumsum(sorted_weights/np.sum(sorted_weights))
      self._confinfo[parname] = {'wfracs': wfracs, 'idxs':idxs, 'Nmask':Nmask}
    idx = np.searchsorted(wfracs, alpha)
    return samples[idxs[idx]]

  def _credible(self,parname,p=None,sigma=None,twoside=False,upper=True):
    # By default, this computes the ETI (equal-tailed-interval)
    if sigma is not None:
      from scipy.special import erf
      if p is not None:
        raise ValueError("Cannot pass both a probability 'p' and a sigma deviation 'sigma'.")
      return self.credible(parname,p=erf(sigma/np.sqrt(2)),twoside=twoside,upper=upper)
    elif p is None:
      raise ValueError("You have to pass either a probability 'p' or a sigma deviation 'sigma'.")
    if p<0 or p>1:
      raise ValueError("Cannot pass {} (p={})".format(("p<0" if p<0 else "p>1"),p))
    if twoside:
      return [self._upper_limit(parname,alpha=(1.-p)/2.),self._upper_limit(parname,alpha=1.-(1.-p)/2.)]
    elif not upper:
      alph = 1-p
    else:
      alph = p
    return self._upper_limit(parname,alpha=alph)
  def credible(self,parnames=None,p=None,sigma=None,twoside=False,upper=True):
    if isinstance(parnames,str):
      return self._credible(parnames,p=p,sigma=sigma,twoside=twoside,upper=upper)
    return {pname:self._credible(pname,p=p,sigma=sigma,twoside=twoside,upper=upper) for pname in (parnames if parnames else self.names[2:])}
  def median(self,parname):
    return self.credible(parname,p=0.5)

  def subdivide(self, subdivisions=None, threshold_min_subdivision=100):
    custom_subdivision = False
    if subdivisions == None:
      subdivisions = 10
    elif not isinstance(subdivisions,int):
      raise Exception("The subdivisions argument for subdivide has to be an integer value. You provided '{}'".format(subdivisions))
    else:
      custom_subdivision = True

    if len(self.lens)>1:
      n_sufficiently_long = np.count_nonzero(self.lens>=threshold_min_subdivision)
      if n_sufficiently_long<1:
        raise Exception("All subdivisions would have sizes below the minimum size allowed '{}'".format(threshold_min_subdivision))
      # If only one of the chains in the folder is long enough, fall back to the single-chain algorithm
      elif n_sufficiently_long == 1:
        custom_subdivision = True
    subfolders = []
    count = 0
    # If only a single folder is present, or if the subdivisions argument has been explicitly set, subdivide differently
    if len(self.lens)==1 or custom_subdivision:
      sublens = np.empty(subdivisions,dtype=int)
      c_per_file = self.N//subdivisions
      remainder = self.N-subdivisions*c_per_file
      if c_per_file < threshold_min_subdivision:
        raise Exception("All subdivisions would have sizes of '{}', while the minimum size allowed is '{}'".format(c_per_file, threshold_min_subdivision))
      for i in range(subdivisions):
        sublens[i] = c_per_file + (1 if remainder>0 else 0)
        if remainder > 0:
          remainder -= 1
      for lvar in sublens:
        subfolders.append(self[count:count+lvar])
        count+=lvar
    else:
      for lvar in self.lens:
        # Ignore small chains with less than threshold_min_subdivision points (we checked above, that there are at least two of these)
        if lvar >= threshold_min_subdivision:
          subfolders.append(self[count:count+lvar])
        count+=lvar
    return subfolders

  def gelman(self, parnames=None, subdivisions=None):
    totN = self.N
    totmean = self.mean(parnames=parnames)
    subfolders = self.subdivide(subdivisions=subdivisions)
    means = np.array([sf.mean(parnames=parnames) for sf in subfolders])
    varias = np.array([np.sqrt(np.diag(sf.cov(parnames=parnames))) for sf in subfolders])
    Ns = np.array([sf.N for sf in subfolders])
    # These should be at least as large as threhsold_min_subdivision of the subdivide method
    W = np.sum(Ns*varias.T,axis=-1)/totN
    B = np.sum(Ns*((means-totmean)**2).T,axis=-1)/(totN-1.)
    R = np.abs(B/W)
    return {n:R[i] for i,n in enumerate((parnames if parnames!=None else self.names[2:]))}

  def max_gelman(self, subdivisions=None):
    subfolders = self.subdivide(subdivisions=subdivisions)
    means = np.array([sf.mean() for sf in subfolders])
    covs = np.array([sf.cov() for sf in subfolders])
    Ns = np.array([sf.N for sf in subfolders])

    mean_of_covs = np.average(covs, weights=Ns, axis=0)
    cov_of_means = np.cov(means.T)#, fweights=Ns) <-- No weighting for increased stability against short outlier chains
    L = np.linalg.cholesky(mean_of_covs)
    Linv = np.linalg.inv(L)
    eigvals = np.linalg.eigvalsh(Linv.dot(cov_of_means).dot(Linv.T))
    return np.max(np.abs(eigvals))

  def to_getdist(self):
    # No logging of warnings temporarily, so getdist won't complain unnecessarily
    #extend using https://github.com/cmbant/getdist/blob/master/getdist/cobaya_interface.py
    import getdist
    getdist.chains.print_load_details = False
    sampler = "mcmc" #We could get this from the log.param file
    names = self.names[2:]
    bounds = self.get_bounds()
    texnames = self._rectify_texnames()
    mcsamples = getdist.MCSamples(
        samples=self.samples.T,
        weights= self.chain._d['N'],
        loglikes=self.chain._d['lnp'],
        names=names,
        labels=texnames,
        sampler=sampler,
        ranges={par:bounds[par] for par in names})
    return mcsamples

  @staticmethod
  def _recursive_rectify(name):
    name = name.replace(" ","")
    idx = name.find("_")
    if idx>=0:
      # Catch only trailing '_' character
      if idx==len(name)-1:
        return name+"{}"
      # Catch almost-trailing "_", always a problem
      if idx==len(name)-2:
        if name[-1]=='_':
          return name[:-2]+r"\_\_"
        elif name[-1]=="{":
          return name[:-2]+"}"
        elif name[-1]=="{":
          return name[:-2]+"{"
        else:
          return name
        return name[:-2]+r"\_"+name[-1]
      # Catch { after underscore, so latex treats everything inside {...} as its own thing
      if name[idx+1]=='{':
        depth = 1
        idxclose = -1
        for i in range(len(name)-idx-2):
          if name[idx+2+i]=='{':
            depth+=1
          elif name[idx+2+i]=='}':
            depth-=1
          if depth==0:
            idxclose = i
            break
        if idxclose>=0:
          subthing = folder._recursive_rectify(name[idx+2:idx+2+idxclose])
          if idx+2+idxclose==len(name)-1:
            return name[:idx+2]+subthing+"}"
          if name[idx+2+idxclose+1]=='_':
            return name[:idx+2]+subthing+r"}\_"+folder._recursive_rectify(name[idx+2+idxclose+2:])
          return name[:idx+2]+subthing+"}"+folder._recursive_rectify(name[idx+2+idxclose+1:])
        else:
          return name[:idx+1]+"{}"+folder._recursive_rectify(name[idx+2:])
      # Catch double underscore -> In this case, they all have to be escaped
      if name[idx+1]=='_':
        return name[:idx]+r"\_\_"+folder._recursive_rectify(name[idx+2:])
      # If none of the above, check if there's another underscore two positions apart
      if name[idx+2]=='_':
        return name[:idx+2]+r"\_"+folder._recursive_rectify(name[idx+3:])
      return name[:idx+1]+folder._recursive_rectify(name[idx+1:])
    else:
      return name

  @staticmethod
  def _rectify_control_characters(name):
    replacedict = {"\n":r"\n","\t":r"\t","\a":r"\a","\r":r"\r"}
    for re, rp in replacedict.items():
      if re in name:
        print(r"liquidcosmo :: /!\ warning, found control character '{}' in your variable name, turning into '{}'. Use r'name' strings to silence this warning!".format(rp,"'+'".join(list(rp))))
        name = name.replace(re,rp)
    return name

  def _rectify_texnames(self):
    return [self._rectify_control_characters(self._recursive_rectify(self._texnames[par])) for par in self.names[2:]]

  def constraint(self,parnames=None):
    if parnames is None:
      parnames = self.names[2:]
    if isinstance(parnames,str):
      return self.constraint([parnames])[parnames]

    # Turn into getdist to get 1d distribution -- TODO :: could also be done via histogramming (if no getdist installed)
    gd = self.to_getdist()
    constraints = {}
    for parname in parnames:
      gd_par = gd.paramNames.parWithName(parname)
      gd_1d = gd.get1DDensity(gd_par.name)
      onesig, twosig = gd_1d.getContourLevels() * self.__limit_safety_factor# or [np.exp(-1**2/2),np.exp(-2**2/2)]

      # Identify whether the posterior sufficiently goes down at the boundary (using getdist)
      lower,upper = self.get_range(parname)
      lower_problem = 0
      upper_problem = 0
      if lower is not None:
        prob = gd_1d(lower)
        # Problem only for 1 sigma
        if prob > onesig:
          lower_problem = 1
        # Problem also for 2 sigma
        elif prob > twosig:
          lower_problem = 2
      if upper is not None:
        prob = gd_1d(upper)
        # Problem only for 1 sigma
        if prob > onesig:
          upper_problem = 1
        # Problem also for 2 sigma
        elif prob > twosig:
          upper_problem = 2

      # With this information, say if there is a constraint
      if lower_problem==1 and upper_problem==1 or (lower_problem==1 and upper_problem==2) or (lower_problem==2 and upper_problem==1):
        #EITHER both a problem at 1 sigma, or one is fine at 1sig and the other cannot give 2 sigma
        constraints[parname] = [[lower,upper],"unconstrained"]
      elif upper_problem==1:
        lowbound =  self._credible(parname,sigma=2,twoside=False,upper=False)
        constraints[parname] = [[lowbound,upper],">"]
      elif lower_problem==1:
        upbound =  self._credible(parname,sigma=2,twoside=False,upper=True)
        constraints[parname] = [[lower,upbound],"<"]
      else:
        lowbound, upbound =  self._credible(parname,sigma=1,twoside=True)
        constraints[parname] = [[lowbound,upbound],"+-"]
    return constraints

  def texconstraints(self,**kwargs):
    return self.texconstraint(**kwargs)
  def texconstraint(self,parnames=None,withdollar=True,withname=True):
    if parnames is None:
      parnames = self.names[2:]
    if isinstance(parnames,str):
      return self.constraint([parnames])[parnames]
    constr = self.constraint(parnames=parnames)
    means = self.mean(parnames=parnames,asdict=True)
    retstr = "\n".join([tex_convert(constr[par],withdollar,means[par],texname=(self._texnames[par] if withname is True else None),name=(par if withname is False else None),equalize=self.__equalize_errors) for par in parnames])
    return retstr

  def plot_getdist(self, **kwargs):
    from .foldercollection import foldercollection
    fc = foldercollection()
    fc.folderlist.append(self)
    return fc.plot_getdist(**kwargs)

  def _add_point(self, spp, add_point, names=None,zorder=None):
    if spp==None or add_point==None:
      return
    if names==None:
      names = self.names
    if isinstance(add_point, dict):
      add_point = [add_point]
    from getdist.plots import ParamInfo
    def __check_arg(arg):
      exception_text = "Argument 'add_point' needs to be a dictionary (or list of dictionaries), containing a value, or containing lists of [value, color] or [mean, sigma, color] for each parameter name -- you provided '{}'.".format(arg)
      if not isinstance(arg,list):
        arg = [arg]
      if len(arg)<1 or len(arg)>3:
        raise Exception(exception_text+" Invalid length!")
      elif len(arg)==1:
        return True, arg[0], 'grey'
      elif len(arg)==2:
        if not isinstance(arg[1],str):
          raise Exception(exception_text+" Invalid color!")
        return True, arg[0], arg[1]
      elif len(arg)==3:
        if not isinstance(arg[2],str):
          raise Exception(exception_text+" Invalid color!")
        return False, arg[0], arg[1], arg[2]
    flag = False
    # Iterate over all plotted subplots
    for i, subplot_arr in enumerate(spp.subplots):
      for j, subplot in enumerate(subplot_arr):
        if subplot != None:
          paramnames_for_subplot = [p.name if isinstance(p, ParamInfo) else p for p in subplot.getdist_params]
          if len(paramnames_for_subplot)>1: # Non-diagonal
            xparam, yparam = paramnames_for_subplot
          else: # Diagonal
            xparam = paramnames_for_subplot[0]
            yparam = None
          for pt in add_point:
            for name in pt:
              arg = __check_arg(pt[name])
              if xparam == name:
                if arg[0]:
                  spp.add_x_marker(arg[1],color=arg[2],ax=[i,j],zorder=zorder)
                else:
                  if yparam:
                    spp.add_x_bands(arg[1],arg[2],color=arg[3],ax=[i,j],zorder=zorder)
                flag = True
              if yparam == name:
                if arg[0]:
                  spp.add_y_marker(arg[1],color=arg[2],ax=[i,j],zorder=zorder)
                else:
                  spp.add_y_bands(arg[1],arg[2],color=arg[3],ax=[i,j],zorder=zorder)
                flag = True
    if not flag:
      raise Exception("Could not find any of '{}' parameters in the generated plot provided as the 'spp' argument.".format([pt.keys() for pt in add_point]))
    return

  def _add_covmat_around_center(self, spp, add_covmat, center_point, names=None,zorder=None,**kwargs):
    if spp is None or add_covmat is None or center_point is None:
      return
    if names is None:
      names = self.names[2:]

    if len(add_covmat)!=len(names):
      raise ValueError("Expected 'add_covmat' of as many entries as names in folder, since no explicit names are provided")
    if len(center_point)!=len(names):
      raise ValueError("Expected 'center_point' of as many entries as names in folder, since no explicit names are provided")

    from getdist.plots import ParamInfo
    flag = False
    # Iterate over all plotted subplots
    for i, subplot_arr in enumerate(spp.subplots):
      for j, subplot in enumerate(subplot_arr):
        if subplot != None:
          paramnames_for_subplot = [p.name if isinstance(p, ParamInfo) else p for p in subplot.getdist_params]
          if len(paramnames_for_subplot)>1: # Non-diagonal
            xparam, yparam = paramnames_for_subplot
            if xparam in names and yparam in names:
              ix, iy = names.index(xparam), names.index(yparam)

              # Create ellipses
              sigma_x = np.sqrt(add_covmat[ix,ix])
              sigma_xy = add_covmat[ix,iy]
              sigma_y = np.sqrt(add_covmat[iy,iy])

              width = 2*np.sqrt(0.5*(sigma_x**2 + sigma_y**2)+np.sqrt((0.5*(sigma_x**2-sigma_y**2))**2+sigma_xy**2))
              height = 2*np.sqrt(0.5*(sigma_x**2 + sigma_y**2)-np.sqrt((0.5*(sigma_x**2-sigma_y**2))**2+sigma_xy**2))
              angle = 0.5*np.arctan(2.*sigma_xy/(sigma_x**2-sigma_y**2)) * 180./np.pi

              ls = kwargs.pop("linestyle","--")
              fill = kwargs.pop("fill",False)
              ecol = kwargs.pop("edgecolor","black")

              from matplotlib.patches import Ellipse
              e1 = Ellipse(xy = [center_point[ix],center_point[iy]], width=width*1.52,height=height*1.52, angle=angle,linestyle=ls,fill=fill,edgecolor=ecol)
              e2 = Ellipse(xy = [center_point[ix],center_point[iy]], width=width*2.48,height=height*2.48, angle=angle,linestyle=ls,fill=fill,edgecolor=ecol)

              subplot.add_patch(e1)
              subplot.add_patch(e2)
              flag=True
    if len(kwargs)>0:
      raise ValueError("Arguments '{}' could not be read.".format(kwargs.keys()))
    if not flag:
      raise KeyError("Could not find any of '{}' parameters in the generated plot provided as the 'spp' argument.".format(add_point.keys()))
    return
  def to_class_dict(self):
    names = self.names[2:]
    samples = self.samples
    if self.logfile == {}:
      arginfo = {}
    else:
      arginfo = self.logfile['arginfo']
    # Define dict before reading so that we can attach things to it during reading of logfile
    retdict = arginfo.copy()
    if self.logfile != {}:
      parinfo = self.logfile['parinfo']
      # Make sure for the class_dict to only copy the cosmo arguments, not nuisance or derived
      newnames= []
      newsamples = []
      for par in parinfo:
        if(parinfo[par]['type']=='cosmo' and parinfo[par]['initialsigma']!=0):
          try:
            ipar = np.nonzero([name==par for name in names])[0][0]
          except IndexError as ie:
            #Very specific fallback for the 100*theta_s -> 100theta_s readability change in MPv3
            if par=='100*theta_s':
              ipar = np.nonzero([name=='100theta_s' for name in names])[0][0]
            else:
              raise
          newnames.append(par)
          newsamples.append(samples[ipar])
        else:
          retdict.update({par:parinfo[par]['initial']})
      names = newnames
      samples = np.array(newsamples).T
    else:
      samples = samples.T
    retlist = []
    retdict.update({name:samples[0][iname] for iname,name in enumerate(names)})
    retlist.append(retdict)
    for i in range(1,self.N):
      retdict.update({name:samples[i][iname] for iname,name in enumerate(names)})
      retlist.append(retdict.copy())
    if self.N==1:
      retlist = retlist[0]
    return retlist

  def to_class_ini(self,filename):
    if self.N != 1:
      raise Exception("Cannot convert a chain file with more than one point to a CLASS ini file")
    inidict = self.to_class_dict()
    with open(filename,"w+") as f:
      for key in inidict:
        f.write("{}={}\n".format(key,inidict[key]))
