import os, re, logging
import use
from ..Platform import platform
from ..Action import Command
from ..Options import Option
from ..File import File
from ..Scanner import CScanner
from ..utils import getarg
from ..conv import to_list, to_iter

class Default(use.Version):
    version = 'default'
    binaries = ['gcc']

    def actions(self, inst, sources, targets=[], options={}):
        return [Command(self.package.options())]

class gcc(use.Package):
    default_binary_filename = 'a.out'
    versions = [Default]

    def __init__(self, *args, **kwargs):
        super(gcc, self).__init__(*args, **kwargs)
        self._opts.binary = 'gcc'
        self._opts.add(Option('compile', '-c'))
        self._opts.add(Option('pic', '-fPIC'))
        self._opts.add(Option('openmp', '-fopenmp'))
        self._opts.add(Option('cxx11', '-std=c++11'))
        self._opts.add(Option('optimise', '-O', space=False))
        self._opts.add(Option('symbols', '-g'))
        if platform.os_name == 'darwin':
            self._opts.add(Option('shared_lib', text='-dynamiclib -install_name {target.abspath} -undefined dynamic_lookup'))
        else:
            self._opts.add(Option('shared_lib', '-shared'))
        self._opts.add(Option('define', '-D', space=False))
        self._opts.add(Option('targets', '-o'))
        self._opts.add(Option('header_dirs', '-I'))
        self._opts.add(Option('library_dirs', '-L'))
        if platform.os_name != 'darwin':
            self._opts.add(Option('rpath_dirs', '-Wl,-rpath=', space=False, abspath=True))
        self._opts.add(Option('sources'))
        self._opts.add(Option('libraries', '-l', space=False))

    ##
    ## gcc's productions. The standard gcc production will
    ## produce a single object file for each source file.
    ##
    def make_productions(self, nodes, inst, opts):
        logging.debug('gcc: Making productions.')

        # Setup the mapping from source extensions to destination extensions.
        if 'suffix' not in opts:
            opts['suffix'] = '.os' if opts.get('pic', False) else '.o'

        # Single target or multitarget?
        single = not opts.get('compile', False)

        # Have the platform order the library directories as appropriate.
        if 'library_dirs' in opts:
            opts['library_dirs'] = platform.order_library_dirs(opts['library_dirs'])

        # Call parent.
        prods = super(gcc, self).make_productions(nodes, inst, opts, single=single)

        # Prepare the scanners for each node.
        for n in nodes:
            if n.scanner is None:
                if os.path.splitext(str(n))[1].lower() in ['.c', '.cc', '.cxx', '.cpp']:
                    n.scanner = CScanner(self.ctx)

        logging.debug('gcc: Done making productions.')
        return prods
