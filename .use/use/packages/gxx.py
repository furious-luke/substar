import gcc

class Default(gcc.Default):
    binaries = ['g++']

class gxx(gcc.gcc):
    versions = [Default]

    def __init__(self, *args, **kwargs):
        super(gxx, self).__init__(*args, **kwargs)
        self.name = 'g++'
        self._opts.binary = 'g++'
