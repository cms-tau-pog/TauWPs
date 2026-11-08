#! /usr/bin/env python3
# Author: Izaak Neutelings (August 2025)
import time
import os
import re
import yaml
from fnmatch import fnmatch
import functools # for @functools.cache

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONFIG_FILE = os.path.join(BASE_DIR, 'config.yaml')
REXP_VAR = re.compile(r"\$\{([^\}]+)\}")

@functools.cache
def loadconfig(cfgname=CONFIG_FILE,verb=0):
  """Load analysis configuration from YAML and normalize."""
  
  # load configuration YAML file
  if verb>=1:
    print(f">>> loadconfig: Opening {cfgname=!r}")
  with open(cfgname, 'r') as infile:
    config = yaml.safe_load(infile) or { }
  
  # sanity checks
  required_keys = [ # required top keys
    'tauAlgos', 'discriminants', 'genmatches',
    'venn', 'dm',
    'title', 'text', 'cut', 'binning',
  ]
  required_subkeys = { # required subkeys of 'text' dict
    'dm': [ 'postfixes', 'postfixes_gm', 'postfixes_wp' ],
    'text': [ 'collection', 'selection', 'algo', 'gm', 'dm', 'process' ]
  }
  if (missing:=[k for k in required_keys if k not in config]):
    raise KeyError(f"loadconfig: Missing keys in {cfgname}: {', '.join(missing)}")
  for key, subkeys in required_subkeys.items():
    if (missing:=[k for k in subkeys if k not in config[key]]):
      raise KeyError(f"loadconfig: Missing keys in {cfgname}:key: {', '.join(missing)}")
  
  # normalize schema
  if verb>=2:
    print(f">>> loadconfig: Normalizing...")
  config['binning'] = {
    k: tuple(v) if isinstance(v, list) else v
    for k, v in config['binning'].items()
  }
  config['text']['dm'] = { int(k): v for k, v in config['text']['dm'].items() }
  if 'GM' not in config['text']:
    config['text']['GM'] = {
      k: (v.capitalize() if (v and v[0] in 'grj') else v)
      for k, v in config['text']['gm'].items()
    }
  
  # expand selection and split strings into lists
  config['text']['selection'] = expand_var_values(config['text']['selection'],verb=verb)
  config['text']['selection'] = {k: v.split('\n') for k, v in config['text']['selection'].items()}
  
  # expand glob patterns of tau algorithms
  keysets = [
    ('dm','postfixes'), ('dm','postfixes_gm'), ('dm','postfixes_wp'),
    ('text','recoAlgo'), ('text','regressAlgo')
  ]
  for key1, key2 in keysets:
    if key2 in config[key1]:
      config[key1][key2] = expand_glob_keys(config[key1][key2],config['tauAlgos'],verb=verb)
    else:
      warning(f"Keys {key1}:{key2} not found! Adding empty...",pre="loadconfig: ")
      config[key1][key2] = { }
  
  # normalize algo labels in text (create algo dict if it doesn't exist)
  for algo in config['tauAlgos']:
    if algo not in config['text']['algo']:
      config['text']['algo'][algo] = algo  # fallback: use algo name as label
  
  # make sure 'dm:postfixes_*' are contained in 'dm:postfixes'
  for key in ['postfixes_gm','postfixes_wp']:
    if key not in config['dm']:
      continue
    for algo in config['dm']['postfixes']:
      if algo not in config['dm'][key]:
        config['dm'][key][algo] = [ ]
        continue
      for postfix in config['dm'][key][algo]:
        if postfix in config['dm']['postfixes'][algo]:
          continue
        warning(f"{postfix=!r} from dm:{key} not in dm:postfixes! Adding...",pre="loadconfig: ")
        config['dm']['postfixes'][algo].append(p)
  
  # add ROC keys: real, elec, muon, jet, fake
  config['rockeys'] = config['genmatches'][:]+['fake'] #+['compare']
  
  return config
  

def expand_glob_keys(old_dict,all_keys,complete=True,verb=0):
  """Expand glob wildcard patterns in given dictionary's keys."""
  new_dict = { }
  unused_keys = all_keys.copy() # avoid overwriticonfig['dm']['gm'].ng
  for old_key, value in old_dict.items():
    for new_key in unused_keys[:]:
      if new_key not in new_dict and fnmatch(new_key,old_key):
        if verb>=2:
          print(f">>> expand_glob_keys: Replacing key {old_key!r} with {new_key!r} for {value=!r}")
        new_dict[new_key] = value
        unused_keys.remove(new_key)
  if complete and unused_keys:
    raise IOError("expand_keys: Not all keys could be expanded!"
                  f" Expanded: {list(new_dict.keys())}, not expanded: {list(unused_keys)}, all keys: {list(all_keys)}, patterns: {list(old_dict.keys())}"
    )
  return new_dict
  

def expand_var_values(old_dict,verb=0):
  """Expand variable of pattern ${varname} in the values of a dictionary."""
  new_dict = old_dict.copy()
  for key in old_dict.keys():
    value = new_dict[key] # use new values in case something changed in the meantime
    for match in REXP_VAR.finditer(value):
      varkey = match.group(1) # variable name
      if varkey in new_dict: # substitute
        if verb>=1:
          print(f">>> expand_var_values: Replacing {match.group(0)!r} in {value!r} with {new_dict[varkey]!r}")
        value = value.replace(match.group(0),new_dict[varkey])
      elif verb>=0:
        warning(f"Could not find {varkey=!r} in {new_dict=}",pre="expand_var_values")
    new_dict[key] = value
  return new_dict
  

def warning(string,pre=''):
  """Print yellow warning."""
  print(f">>> \033[93m{pre}\033[1mWARNING!\033[0m\033[93m {string}\033[0m")
  

def elapsed(time0):
  """Format elapsed time as hours, minutes, seconds."""
  secs = time.perf_counter()-time0
  mins, secs = divmod(secs,60)
  hours, mins = divmod(mins,60)
  if hours: return f"{hours:.0f}h {mins:.0f}m {secs:.0f}s"
  elif mins: return f"{mins:.0f}m {secs:.2g}s"
  return f"{secs:.3g}s"
  
