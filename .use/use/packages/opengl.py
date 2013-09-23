import use

class Default(use.Version):
    version = 'default'
    headers = ['GL/gl.h']
    libraries = ['GL']

class opengl(use.Package):
    versions = [Default]
