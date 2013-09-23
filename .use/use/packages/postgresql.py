import use

class Default(use.Version):
    version = 'default'
    headers = ['libpq-fe.h']
    libraries = ['pq']

class postgresql(use.Package):
    versions = [Default]
