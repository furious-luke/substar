import use

class postgresql(use.Feature):
    dependencies = ['postgresql']
    headers = ['soci/postgresql/soci-postgresql.h']
    libraries = ['soci_postgresql']

class sqlite3(use.Feature):
    dependencies = ['sqlite3']
    headers = ['soci/sqlite3/soci-sqlite3.h']
    libraries = ['soci_sqlite3']

class Default(use.Version):
    version = 'default'
    headers = ['soci/soci.h']
    libraries = ['soci_core']
    features = [postgresql, sqlite3]

class soci(use.Package):
    versions = [Default]
    header_sub_dirs = ['include/soci']
    # optional_dependencies = ['cmake']
    url = 'http://downloads.sourceforge.net/project/soci/soci/soci-3.1.0/soci-3.1.0.zip'

    def build_handler(self):

        # SOCI can use Boost, so check to see if Boost is in our
        # set of configuration options and set accordingly.
        cmake = 'cmake -DCMAKE_INSTALL_PREFIX:PATH={prefix}'
        # boost = config.package(config.packages.boost)
        # if boost and boost.found and boost.base_dir:
        #     cmake += ' -DBOOST_DIR:PATH=' + boost.base_dir
        cmake += ' -DWITH_BOOST=off'

        # Turn on release mode to prevent a very odd bug
        # in SOCI that produces an assertion whenever you
        # try to read data from PostgreSQL into a vector
        # larger than numeric_limits<unsigned short>::max().
        cmake += ' -DCMAKE_BUILD_TYPE:STRING=release'

        # Check for sqlite3, like boost.
        sqlite = config.package(config.packages.sqlite3)
        if sqlite and sqlite.found and sqlite.base_dir:
            cmake += ' -DSQLITE3_INCLUDE_DIR:PATH=' + sqlite.include_directories()
            cmake += ' -DSQLITE3_LIBRARY:FILEPATH=' + sqlite.libraries()

        # Check for MySQL.
        pkg = config.package(config.packages.MySQL)
        if pkg and pkg.found and pkg.base_dir:
            cmake += ' -DMYSQL_INCLUDE_DIR:PATH=' + pkg.include_directories()
            cmake += ' -DMYSQL_LIBRARY:FILEPATH=' + pkg.libraries()

        # Check for PostgreSQL.
        pkg = config.package(config.packages.PostgreSQL)
        if pkg and pkg.found and pkg.base_dir:
            cmake += ' -DPOSTGRESQL_INCLUDE_DIR:PATH=' + pkg.include_directories()
            cmake += ' -DPOSTGRESQL_LIBRARY:FILEPATH=' + pkg.libraries()

        # For some reason SOCI is incompatible with gcc 4.7.1.
        # Need to switch off testing and the empty thingy.
        cmake += ' -DSOCI_TEST:BOOL=off -DSOCI_EMPTY:BOOL=off'

        cmake += ' .'
        return [cmake, 'make', 'make install']
