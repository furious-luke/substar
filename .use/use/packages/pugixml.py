import use

class Default(use.Version):
    version = 'default'
    headers = ['pugixml.hpp']
    libraries = ['pugixml']

class pugixml(use.Package):
    default_target_node = use.File
    default_builder = use.Builder
    versions = [Default]
    url = 'http://pugixml.googlecode.com/files/pugixml-1.2.tar.gz'

    def build_handler(self):
        return ['cmake scripts -DBUILD_SHARED_LIBS=yes -DCMAKE_INSTALL_PREFIX={prefix}',
                'make',
                'make install']
