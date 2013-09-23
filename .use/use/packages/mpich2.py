import os, re, logging
import use
from ..utils import getarg
from ..conv import to_list

##
## A feature for mpicc.
##
class mpicc(use.Feature):
    binaries = ['mpicc']

    def __init__(self, *args, **kwargs):
        super(mpicc, self).__init__(*args, **kwargs)

        # We will need to figure out which compiler type this mimics.
        self.cc = None

    def actions(self, sources, targets=[]):
        return self.cc.expand(*args, **kwargs)

    def expand(self, *args, **kwargs):
        return self.cc.expand(*args, **kwargs)

class pmpich(use.Feature):
    libraries = ['pmpich']

##
##
##
class Default(use.Version):
    version = 'default'
    headers = ['mpi.h']
    libraries = ['mpich']
    features = [pmpich, mpicc]
    url = 'http://www.mcs.anl.gov/research/projects/mpich2staging/goodell/downloads/tarballs/mpich2-current.tar.gz'

##
## MPICH2
##
class mpich2(use.Package):
    default_target_node = use.File
    default_builder = use.Builder
    versions = [Default]
    environ_name = 'MPI'

    def build_handler(self):
        return ['./configure --prefix={prefix} --enable-shared --disable-fc --disable-f77',
                'make',
                'make install']
