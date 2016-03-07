from __future__ import division

# DIALS version numbers are constructed from
#  1. a common prefix
__dials_version_format = "DIALS %s"
#  2. the most recent annotated git tag (or failing that: a default string)
__dials_version_default = "1.dev"
#  3. a dash followed by the number of commits since that tag
#  4. a dash followed by a lowercase 'g' and the current commit id

# When run from a development installation the version information is extracted
# from the git repository. Otherwise it is read from the file '.gitversion' in the
# DIALS module directory.

def dials_version():
  '''Try to obtain the current git revision number
     and store a copy in .gitversion'''
  version = None

  try:
    import libtbx.load_env
    import os
    dials_path = libtbx.env.dist_path('dials')
    version_file = os.path.join(dials_path, '.gitversion')

    # 1. Try to access information in .git directory
    #    Regenerate .gitversion if possible
    if os.path.exists(os.path.join(dials_path, '.git')):
      try:
        import subprocess
        with open(os.devnull, 'w') as devnull:
          version = subprocess.check_output(
            ["git", "describe", "--long"], cwd=dials_path, stderr=devnull).rstrip()
          if version[0] == 'v':
            version = version[1:].replace('.0-','-')
          version = version.replace('-', '.', 1)
          try:
            branch = subprocess.check_output(["git", "describe", "--contains", "--all", "HEAD"], cwd=dials_path, stderr=devnull).rstrip()
            if 'release' in branch:
              version = version + '-release'
          except Exception:
            pass
        with open(version_file, 'w') as gv:
          gv.write(version)
      except Exception:
        if version == "": version = None

    # 2. If .git directory missing or 'git describe' failed, read .gitversion
    if (version is None) and os.path.exists(version_file):
      with open(version_file, 'r') as gv:
        version = gv.read().rstrip()
  except Exception:
    pass

  if version is None:
    version = __dials_version_format % __dials_version_default
  else:
    version = __dials_version_format % version

  return version
