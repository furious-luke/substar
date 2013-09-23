import use
from ..Platform import platform

# class regex(use.Feature):
#     libraries = ['boost_regex']

# class iostreams(use.Feature):
#     libraries = ['boost_iostreams']

# class filesystem(use.Feature):
#     libraries = ['boost_filesystem']

class Default(use.Version):
    version = 'default'
    headers = ['boost/optional.hpp']
    libraries = [['boost_system', 'boost_regex', 'boost_iostreams', 'boost_filesystem'],
                 ['boost_system-mt', 'boost_regex-mt', 'boost_iostreams-mt', 'boost_filesystem-mt']]
    # features = [regex, iostreams, filesystem]

class boost(use.Package):
    default_builder = use.Builder
    versions = [Default]
