import os
from .conv import to_list

class Location(object):

    def __init__(self, base, bin_dirs, hdr_dirs, lib_dirs):
        self.base = base
        self.binary_dirs = list(self._iter_dirs(bin_dirs)) if bin_dirs is not None else []
        self.header_dirs = list(self._iter_dirs(hdr_dirs)) if hdr_dirs is not None else []
        self.library_dirs = list(self._iter_dirs(lib_dirs)) if lib_dirs is not None else []

    def __eq__(self, op):
        if self.base != op.base:
            return False
        for dirs in ['binary_dirs', 'header_dirs', 'library_dirs']:
            for d in getattr(self, dirs):
                if d not in getattr(op, dirs):
                    return False
        return True

    def __repr__(self):
        return self.text().replace('\n', ', ')

    def text(self, indent=0):
        text = []
        if self.base:
            text.append('base: ' + self.base)
        if self.binary_dirs:
            text.append('binary dirs: ' + str(self.binary_dirs))
        if self.header_dirs:
            text.append('header dirs: ' + str(self.header_dirs))
        if self.library_dirs:
            text.append('library dirs: ' + str(self.library_dirs))
        return ('\n' + indent*' ').join(text)

    def _iter_dirs(self, dirs):
        dirs = to_list(dirs)
        for d in dirs:
            if os.path.isabs(d):
                path = d
            else:
                path = os.path.join(self.base, d)
            if os.path.exists(path):
                yield path
