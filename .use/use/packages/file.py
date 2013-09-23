import re
import use
from ..utils import getarg
from ..conv import to_list

##
##
##
class Builder(use.Builder):

    _prog = re.compile('(' + ')|('.join([
        'symbolic link to',
        'ar archive',
        'shared object',
        'executable,',
        'C program text',
        r'C\+\+ program text',
        'Python script',
    ]) + ')')

    def _parse(self):
        assert len(self.sources) == 1 and len(self.targets) == 1, 'Can\'t run "file" on multiple sources/targets.'
        match = self._prog.search(self.actions[0].stdout)
        assert match and match.lastindex > 0
        self.targets[0].file_type = match.lastindex - 1

##
##
##
class FileType(use.Node):

    SYMBOLIC_LINK = 0
    AR_ARCHIVE = 1
    SHARED_OBJECT = 2
    EXECUTABLE = 3
    ASCII_C_TEXT = 4
    ASCII_CXX_TEXT = 5
    ASCII_PYTHON_TEXT = 6

    def __init__(self):
        super(FileType, self).__init__()
        self.file_type = None

##
##
##
class Default(use.Version):
    version = 'default'
    binaries = ['file']

    def _actions(self, sources, targets=[]):
        return [use.Action('%s %s'%(self.binaries()[0].path, sources[0].path))]

##
## GNU "file" tool.
##
class file(use.Package):
    default_target_node = FileType
    default_builder = Builder
    versions = [Default]
