import os, platform, sys, glob
from .Location import Location
from .conv import to_list
from .utils import strip_missing, run_command

class Platform(object):

    def __init__(self):
        self.base_dirs = []
        self.search_dirs = []
        self.binary_sub_dirs = []
        self.header_sub_dirs = []
        self.library_sub_dirs = []
        self.system_header_dirs = []
        self.system_library_dirs = []
        self._static_lib_str = 'lib{name}.a'
        self._shared_lib_str = 'lib{name}.so'
        self._set_os(platform.system())

    def iter_locations(self, patterns=[], hdr_sub_dirs=[]):
        hdr_sub_dirs = [hdr_sub_dirs + d for d in self.header_sub_dirs]
        patterns = to_list(patterns)
        for base in self.base_dirs:
            for loc in self.iter_base_locations(base, None, hdr_sub_dirs):
                yield loc
        for sd in self.search_dirs:
            for pattern in patterns:
                for base in glob.iglob(os.path.join(sd, pattern)):
                    for loc in self.iter_base_locations(base, None, hdr_sub_dirs):
                        yield loc

    def iter_base_locations(self, base_dir, bin_dir=None, hdr_dir=None, lib_dir=None):
        bin_sub_dirs = to_list(bin_dir) if bin_dir is not None else self.binary_sub_dirs
        hdr_sub_dirs = to_list(hdr_dir) if hdr_dir is not None else self.header_sub_dirs
        lib_sub_dirs = to_list(lib_dir) if lib_dir is not None else self.library_sub_dirs
        for bin_dirs in bin_sub_dirs:
            for hdr_dirs in hdr_sub_dirs:
                for lib_dirs in lib_sub_dirs:
                    yield Location(base_dir, bin_dirs, hdr_dirs, lib_dirs)

    def make_static_library(self, name):
        return os.path.join(os.path.dirname(name), self._static_lib_str.format(name=os.path.basename(name)))

    def make_shared_library(self, name):
        return os.path.join(os.path.dirname(name), self._shared_lib_str.format(name=os.path.basename(name)))

    def order_library_dirs(self, lib_dirs):
        dirs = set(['/usr/lib/x86_64-linux-gnu', '/usr/lib64/x86_64-linux-gnu',
                    '/usr/lib/i386-linux-gnu', '/usr/lib64/i386-linux-gnu'])
        new_lib_dirs = []
        post = []
        for d in lib_dirs:
            if d in dirs:
                post.append(d)
            else:
                new_lib_dirs.append(d)
        new_lib_dirs.extend(post)
        return new_lib_dirs

    def _set_os(self, os_name):
        self.os_name = os_name.lower()
        try:
            getattr(self, '_' + self.os_name)()
        except AttributeError:
            print 'Uh oh, the platform you\'re running on is not yet supported.'
            print 'Please contact "furious.luke@gmail.com" to fix this.'
            sys.exit(1)

    def _darwin(self):
        self._linux()
        self._shared_lib_str = 'lib{name}.dylib'

    def _linux(self):
        self._unix()

    def _unix(self):
        mac = '-'.join([platform.machine(), 'linux', 'gnu'])
        self.base_dirs = strip_missing(['/usr', '/opt', '/opt/local'])
        self.search_dirs = strip_missing(['/usr/local', '/opt/local', os.environ['HOME'], os.path.join(os.environ['HOME'], 'soft')])
        self.binary_sub_dirs = [['bin']]
        self.header_sub_dirs = [['include']]
        self.library_sub_dirs = [['lib', 'lib/' + mac], ['lib64', 'lib64/' + mac]]
        self.system_header_dirs = ['/usr/include']
        self.system_library_dirs = ['/usr/lib', '/usr/lib64']

platform = Platform()
