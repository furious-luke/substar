import use

class default(use.Version):
    headers = ['regexp/regexp.hh']
    libraries = ['regexp']

class regexp(use.Package):
    versions = [default]
    dependencies = ['libhpc']
