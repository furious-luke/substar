import use

class mpicxx(use.Feature):
    name = 'mpic++'
    binaries = ['mpic++']
    libraries = ['mpi_cxx']

class Default(use.Version):
    version = 'default'
    headers = ['mpi.h']
    libraries = ['mpi']
    features = [mpicxx]

class openmpi(use.Package):
    versions = [Default]
    environ_name = 'MPI'
