import use

class Default(use.Version):
    version = 'default'

class identity(use.Package):
    versions = [Default]

    def expand(self, *args, **kwargs):
        return None

    def apply(self, *args, **kwargs):
        pass
