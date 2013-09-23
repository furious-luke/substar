import use

class default(use.Version):
    headers = ['openssl/ssl.h']
    libraries = ['ssl']

class openssl(use.Package):
    versions = [default]
