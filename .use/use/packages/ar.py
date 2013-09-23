import use
from ..Action import Command
from ..Options import Option

class Default(use.Version):
    version = 'default'
    binaries = ['ar']

    def actions(self, inst, sources, targets=[], options={}):
        return [Command(self.package.options())]

class ar(use.Package):
    default_target_node = use.File
    default_builder = use.Builder
    versions = [Default]

    def __init__(self, *args, **kwargs):
        super(ar, self).__init__(*args, **kwargs)
        self._opts.binary = 'ar'
        self._opts.add(Option('add', 'rs'))
        self._opts.add(Option('targets'))
        self._opts.add(Option('sources'))

    def make_productions(self, nodes, inst, opts):
        return super(ar, self).make_productions(nodes, inst, opts, single=True)
