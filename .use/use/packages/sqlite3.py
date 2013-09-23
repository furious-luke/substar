import use

class default(use.Version):
    headers = ['sqlite3.h']
    libraries = ['sqlite3']

class sqlite3(use.Package):
    versions = [default]
