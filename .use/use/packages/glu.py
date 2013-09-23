import use

class Default(use.Version):
    version = 'default'
    headers = ['GL/glu.h']
    libraries = ['GLU']

class glu(use.Package):
    versions = [Default]
    dependencies = ['opengl']
