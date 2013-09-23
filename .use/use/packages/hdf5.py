import os
import use

class parallel(use.Feature):
    dependencies = ['mpi']

    def check(self, inst):
        for lib_dir in inst.library_dirs:
            with open(os.path.join(lib_dir, 'libhdf5.settings'), 'r') as setf:
                data = setf.read()
                if data.find('Parallel HDF5: yes') != -1:
                    return True
        return False

class Default(use.Version):
    version = 'default'
    headers = ['hdf5.h']
    libraries = ['hdf5']
    features = [parallel]

class hdf5(use.Package):
    default_target_node = use.File
    default_builder = use.Builder
    versions = [Default]
