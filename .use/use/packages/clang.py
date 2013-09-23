import gcc

class Default(gcc.Default):
    binaries = ['clang']

class clang(gcc.gcc):
    versions = [Default]

    def __init__(self, *args, **kwargs):
        super(clang, self).__init__(*args, **kwargs)
        self._opts.binary = 'clang'
