import use
from ..Platform import platform

class cxx_compiler(use.Package):
    name = 'C++ compiler'
    option_name = 'cxx'
    sub_packages = ['gxx', 'clangxx']

    def iter_sub_packages(self):

        # If we're building on Darwin, favor clang ahead of
        # GNU gcc or any other gcc.
        if platform.os_name == 'darwin':
            yield self._sub_pkg_map['clangxx']
            for n in cxx_compiler.sub_packages:
                if n != 'clangxx':
                    yield self._sub_pkg_map[n]

        else:
            for n in super(cxx_compiler, self).iter_sub_packages():
                yield n
