import numpy as np
_round_exponential_max = 5
_round_exponential_min = -5
def round_reasonable(val, errp=None,errm=None , digits=1):
  # Do some accounting a priori
  if digits<0:
    raise ValueError("Cannot print negative number of significant digits")
  if errm is not None and errp is None:
    raise ValueError("Cannot give only negative but no positive error bar")
  # Check for arrays instead of values, convert individually
  if isinstance(val,(list,tuple,np.ndarray)):
    if errp and not isinstance(errp,(list,tuple,np.ndarray)) and not len(errp)==len(val):
      raise ValueError("Cannot pass a list of vals but only a single errp")
    if errm and not isinstance(errm,(list,tuple,np.ndarray)) and not len(errp)==len(val):
      raise ValueError("Cannot pass a list of vals but only a single errm")
    return [round_reasonable(v,errp=(errp[iv] if errp else None),errm=(errm[iv] if errm else None),digits=digits) for iv,v in enumerate(val)]
  # Only value, no errors
  if errp is None:
    if val==0:
      first_sigfig = 0
    else:
      first_sigfig = int(np.floor(np.log10(np.abs(val))))
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      return "{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \cdot 10^{"+str(first_sigfig)+"}"
    else:
      return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)
  # Symmetric error
  elif errm is None:
    if errp<0:
      raise ValueError("Cannot pass negative errp")
    # error = 0 is equivalent to no error
    if errp==0:
      return round_reasonable(val, errp=None,errm=None , digits=digits)
    # Use error sigfig to do rest
    first_sigfig = int(np.floor(np.log10(np.abs(errp))))
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      return "("+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \pm "+"{:.{acc}f}".format(errp*10**(-first_sigfig),acc=digits)+") \cdot 10^{"+str(first_sigfig)+"}"
    else:
      return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)+" \pm "+"{:.{acc}f}".format(errp,acc=-first_sigfig+digits)
  else:
    if errp<0:
      raise ValueError("Cannot pass negative errp")
    if errm<0:
      errm=-errm
    if errp==0 and errm==0:
      return round_reasonable(val, errp=None,errm=None , digits=digits)
    if errp==0:
      errp_sigfig = 0
    else:
      errp_sigfig = int(np.floor(np.log10(np.abs(errp))))
    if errm==0:
      errm_sigfig = 0
    else:
      errm_sigfig = int(np.floor(np.log10(np.abs(errm))))
    # Choose the LARGER error for the significant digit
    first_sigfig = max(errp_sigfig,errm_sigfig)
    if first_sigfig < _round_exponential_min or first_sigfig > _round_exponential_max:
      pluserrstr = "{:.{acc}f}".format(errp*10**(-first_sigfig),acc=digits)
      minuserrstr = "{:.{acc}f}".format(errm*10**(-first_sigfig),acc=digits)
      if pluserrstr!=minuserrstr:
        return "({"+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+"}^{+"+pluserrstr+"}_{-"+minuserrstr+"}"+") \cdot 10^{"+str(first_sigfig)+"}"
      else:
        return "("+"{:.{acc}f}".format(val*10**(-first_sigfig),acc=digits)+" \pm "+pluserrstr+") \cdot 10^{"+str(first_sigfig)+"}"
    else:
      pluserrstr = "{:.{acc}f}".format(errp,acc=-first_sigfig+digits)
      minuserrstr = "{:.{acc}f}".format(errm,acc=-first_sigfig+digits)
      if pluserrstr!=minuserrstr:
        return "{"+"{:.{acc}f}".format(val,acc=-first_sigfig+digits)+"}^{+"+pluserrstr+"}_{-"+minuserrstr+"}"
      else:
        return "{:.{acc}f}".format(val,acc=-first_sigfig+digits)+" \pm "+"{:.{acc}f}".format(errp,acc=-first_sigfig+digits)
