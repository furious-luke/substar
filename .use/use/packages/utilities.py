import use

class default(use.Version):
    headers = ['utilities/utilities.hh']
    libraries = ['utilities']

class utilities(use.Package):
    versions = [default]
    dependencies = ['libhpc']
