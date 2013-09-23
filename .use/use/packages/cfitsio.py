import use

class Default(use.Version):
    version = 'default'
    libraries = ['cfitsio']

class cfitsio(use.Package):
    versions = [Default]
