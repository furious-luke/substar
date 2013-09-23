import use

class Default(use.Version):
    version = 'default'
    headers = ['GL/glut.h']
    libraries = ['glut']

class glut(use.Package):
    versions = [Default]
    dependencies = ['glu']
