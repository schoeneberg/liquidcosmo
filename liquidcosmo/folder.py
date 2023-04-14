import os
import numpy as np
import signal
import multiprocessing
from .chain import chain
import ast
from datetime import date
from collections import OrderedDict
from copy import deepcopy
from .matplotlib_defaults import default_settings

class _lq_code_type:
  montepython = 0
  cobaya = 1

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

#from collections import OrderedDict
class folder:
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
    self.path = None

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
    a.path = self.path
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
    a.path = deepcopy(self.path)
    return a

  # Construct a folder object (containing multiple physical chains) from a path
  @classmethod
  def load(obj, path, kind="all", burnin_threshold = 3, verbose = 0):
    a = obj()
    a.verbose = verbose
    a._foldername, a._allchains, a._chainprefix, a._code = obj._resolve_chainname(path, kind=kind)
    a.lens = None
    a._texnames = None
    a._arr = None
    a._narr = None
    a._log = None
    a.path = os.path.abspath(a._foldername)
    a.get_chain(burnin_threshold=burnin_threshold)
    return a

  # -- Load a given chain (for a given full filename)
  def __ldchain__(filename,verbose=False):
    if verbose:
      print("liquidcosmo :: Loading chain {}".format(filename))
    arr = np.genfromtxt(filename).T
    return arr

  # -- Resolve the chain names of chains --
  @classmethod
  def _resolve_chainname(obj,path,kind="all"):
    if os.path.isfile(path):
      folder = os.path.dirname(path)
      allchains = [path]
    elif os.path.isdir(os.path.join("chains",path)):
      folder = os.path.join("chains",path)
      allchains, code = obj._get_chainnames(folder,kind=kind)
    elif os.path.isdir(path):
      folder = path
      allchains, code = obj._get_chainnames(folder,kind=kind)
      prefix = obj._resolve_prefix(allchains[0])
    elif os.path.exists(os.path.join("chains",path)):
      folder = os.path.dirname(os.path.join("chains",path))
      allchains, code = obj._get_chainnames(folder,kind=kind)
      prefix = obj._resolve_prefix(allchains[0])
    elif os.path.exists(os.path.join("chains",path+"__1.txt")):
      folder = os.path.dirname(os.path.join("chains",path+"__1.txt"))
      allchains, code = obj._get_chainnames(folder,kind=kind)
      allchains = [ch for ch in allchains if path in ch]
      prefix = obj._resolve_prefix(allchains[0])
    else:
      folder = os.path.dirname(path)
      if folder=="":
        raise Exception("Could not find a chain from : "+path)
      allchains = obj._get_chainnames(folder,kind=kind)
    prefix = obj._resolve_prefix(allchains[0])
    return folder,allchains,prefix,code

  @classmethod
  def _resolve_prefix(obj, path):
    if "__" in path:
      return path.split("__")[0]
    elif ".txt" in path:
      return "".join(path.split(".")[:-2])
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
    else:
      code = _lq_code_type.cobaya
      chain_targets = [chain for chain in allfilenames if ".txt" in chain]

    if kind=="newest" and code== _lq_code_type.montepython:
      chain_targets = np.sort(chain_targets)[::-1]
      newest_chain_dates = chain_targets[0].split("__")[0]
      allchains= [os.path.join(directory,chain) for chain in chain_targets if newest_chain_dates in chain]

    elif kind=="all":
      allchains= [os.path.join(directory,chain) for chain in chain_targets]

    else:
      raise Exception("Kind can either be 'all' or 'newest'")

    return allchains, code

  # -- Load all points from folder --
  def _get_array(self,excludesmall=True,burnin_threshold=3):
    if self._arr is not None:
      return self._arr

    chainnames = self._allchains

    if excludesmall:
      chainnames = [chain for chain in chainnames if not "_1__" in chain]

    sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)
    pool = multiprocessing.Pool(8)
    signal.signal(signal.SIGINT, sigint_handler)
    try:
      filearr = pool.map_async(folder.__ldchain__, chainnames)
      filearr = filearr.get(60) # 60 seconds (timeout)
    except KeyboardInterrupt:
      print("You stopped the loading of the chains by KeyboardInterrupt")
      pool.terminate()
    else:
      pool.close()
    pool.join()

    filearr = [fa for fa in filearr if fa!=[] and fa.ndim>1]
    if(len(filearr))==0:
      raise Exception("There is probably a problem with the chain folder '{}' that is attempted to be analyzed. Please make sure that the folder is non-empty and the chains are not empty files.".format(self._foldername))
    arrs = [[] for i in range(len(filearr[0]))]

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


    for iparam in range(len(filearr[0])):
      for j in range(len(filearr)):
        arrs[iparam] = np.concatenate([arrs[iparam],filearr[j][iparam][:]])
    arrs = np.array(arrs)

    self._arr = arrs
    return arrs

  # -- Get dictionary of all parameters in .paramnames file --
  # Returns a dictionary with parameter name and index
  def _load_names(self):
    retdict = {}
    texdict = {'N':'Multiplicity','lnp':'-\\ln(\\mathcal{L})'}
    index = 2
    if self._code == _lq_code_type.montepython:
      with open(self._chainprefix+"_.paramnames") as parnames:
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
        names = [x for x in line[1:].split(" ") if x]
        retdict = {name:i for i,name in enumerate(names)}
        texdict = {name:name for name in names}
    else:
      raise Exception("Missing code type (report to developer)")
    return retdict,texdict

  # Create the chain object (equivalent to a named dictionary or arrays)
  def get_chain(self,excludesmall=True,burnin_threshold=5):
    if not (self._narr is None):
      return self._narr
    arr = self._get_array(excludesmall=excludesmall,burnin_threshold=burnin_threshold)
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
      print(e)
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
          par['min'] = par['prior'].pop('loc',0)
          par['max'] = par['min']+par['prior'].pop('scale',1)
        if 'min' in par:
          par['bound'][0]  = par['min']
        if 'max' in par:
          par['bound'][1]  = par['max']
      arginfo = logdict.pop('theory')
      loginfo = logdict
      self._log = {'loginfo':loginfo,'parinfo':parinfo,'arginfo':arginfo,'lklopts':lklopts}
    return self._log

  def write(self, fname, codetype="montepython"):
    if not (codetype=="montepython" or codetype=="Montepython" or codetype=="MontePython" or codetype=="MONTEPYTHON"):
      raise Exception("Can currently only write in MontePython style")
    if self._code != _lq_code_type.montepython:
      raise Exception("Sadly this conversion from code type '{}' to 'montepython' is not yet supported".format(self._code))
    if os.path.exists(fname):
      raise Exception("File path already exists : ",fname)
    else:
      # How many points to put into the individual chains?
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
      os.mkdir(fname)
      index = 0
      for i, c in enumerate(counts):
        with open(os.path.join(fname,str(date.today())+"_"+str(c)+"__"+str(i)+".txt"),"w") as ofile:
          arr = self.chain[index:index+c+1]._d
          index+=c
          ofile.write("# "+"\t".join(arr.keys())+"\n")
          for j in range(len(arr['N'])):
            ofile.write("%.4g"%arr['N'][j])
            ofile.write(" %.6g\t"%arr['lnp'][j])
            for k in list(arr.keys())[2:]:
              ofile.write("%.6e\t"%arr[k][j])
            ofile.write("\n")
      with open(os.path.join(fname,str(date.today())+"_"+str(c)+"_.paramnames"),"w") as ofile:
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

  def set_range(self,parname,lower=None,upper=None):
    if self.logfile != {} and parname in self.logfile['parinfo']:
      if isinstance(lower,list) or type(lower) is np.ndarray:
        if upper:
          raise Exception("Cannot have 'upper!=None' but also 'lower' be a list")
        elif len(lower)!=2:
          raise Exception("If you provide a list, it must be exactly 2 elements long")
        else:
          upper = lower[1]
          lower = lower[0]
      self.logfile['parinfo'][parname]['bound'] = [lower,upper]
  def get_range(self,parname):
    if self.logfile != {} and parname in self.logfile['parinfo']:
      return self.logfile['parinfo'][parname]['bound']
    else:
      return [None,None]
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

  def mean(self,parnames=None):
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

  def _parmean(self,parname):
    return np.average(self[parname],weights=self['N'])
  def _parcov(self,parnames):
    if self.N==1:
      raise Exception("Cannot compute covariance of a chain with only a single point!")
    return np.cov([self[parname] for parname in parnames],fweights=self['N'])

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
    mcsamples = getdist.MCSamples(
        samples=self.samples.T,
        weights= self.chain._d['N'],
        loglikes=self.chain._d['lnp'],
        names=names,
        labels=[self._texnames[par] for par in names],
        sampler=sampler,
        ranges={par:bounds[par] for par in names})
    return mcsamples

  def plot_getdist(self, ax=None,color=None,add_point=None,**kwargs):
    from getdist.plots import get_subplot_plotter
    gdfolder = self.to_getdist()
    spp = get_subplot_plotter(settings=default_settings)
    if 'filled' not in kwargs:
      kwargs['filled'] = True
    spp.triangle_plot([gdfolder],colors=[color],
      line_args=({'color':c} if color else None), **kwargs)
    self._add_point(spp,add_point)

  def _add_point(self, spp, add_point, names=None,zorder=None):
    if spp==None or add_point==None:
      return
    if names==None:
      names = self.names
    from getdist.plots import ParamInfo
    def __check_arg(arg):
      exception_text = "Argument 'add_point' needs to be a dictionary, containing a value, or containing lists of [value, color] or [mean, sigma, color] -- you provided '{}'.".format(arg)
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
          for name in add_point:
            arg = __check_arg(add_point[name])
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
      raise Exception("Could not find any of '{}' parameters in the generated plot provided as the 'spp' argument.".format(add_point.keys()))
    return

  def to_class_dict(self):
    names = self.names[2:]
    samples = self.samples
    if self.logfile != {}:
      arginfo = self.logfile['arginfo']
      parinfo = self.logfile['parinfo']
      # Make sure for the class_dict to only copy the cosmo arguments, not nuisance or derived
      newnames= []
      newsamples = []
      for par in parinfo:
        if(parinfo[par]['type']=='cosmo'):
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
      names = newnames
      samples = np.array(newsamples).T
    else:
      arginfo = {}
      samples = samples.T
    retlist = []
    # By defining the dict here, we can reuse it in the loop, saving computational time
    retdict = arginfo.copy()
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
