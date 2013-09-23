import sys, os, argparse, pickle
from utils import getarg, load_class
from .Use import Use
from .Rule import *
from .Resolver import Resolver
from .Argument import Arguments
from .Options import OptionDict
from .File import File
from .conv import to_list
import logging

##
## Wraps a single build context. The context manages the
## overall run, including commencement and storage. (??)
##
class Context(object):

    def __init__(self):
        self.packages = []
        self.rules = []
        self.uses = []
        self._node_map = {}
        self.targets = []
        self.resolver = Resolver()
        self._pkg_map = {}
        self.crcs = {}
        self.src_crcs = {}
        self.old_bldrs = {}

        self.arguments = None
        self._def_args = {}
        self._arg_map = {}
        self.parser = argparse.ArgumentParser('"Use": Software configuration and build.')
        self.parser.add_argument('targets', nargs='*', help='Specify build targets.')
        self.parser.add_argument('-s', dest='show_config', action='store_true', help='Show current configuration.')
        self.new_arguments()('--enable-download-all', dest='download_all', action='boolean', help='Download and install all dependencies.')

    def __eq__(self, op):

        # Must have the same packages.
        for pkg in self.packages:
            if pkg not in op.packages:
                return False

        # # Must have same rules.
        # for rule in self.rules:
        #     if rule not in op.rules:
        #         return False

        # Must have same uses.
        for use in self.uses:
            if use not in op.uses:
                return False

        # # Must be using the same resolver.
        # if self.resolver != op.resolver:
        #     return False

        return True

    def __ne__(self, op):
        return not self.__eq__(op)

    ##
    ##
    ##
    def parse_arguments(self):

        # Add arguments.
        for pkg in self.packages:
            pkg.add_arguments(self.parser)

        # Parse.
        self.arguments = self.parser.parse_args()

        # Check if we have an old structure to use, unless the user
        # requested a reconfiguration.
        if 'configure' not in self.arguments.targets and os.path.exists('.use.db'):
            with open('.use.db', 'r') as inf:
                old_ctx = pickle.load(inf)

            # Only update arguments if allowed.
            for k, v in self.arguments.__dict__.iteritems():
                if v is not None:
                    old_ctx.arguments.__dict__[k] = v
            self.arguments = old_ctx.arguments

        # Parse the options now.
        for use in self.uses:
            if use.options is not None:
                use.options.parse(self)
        for rule in self.rules:
            if rule.options is not None:
                rule.options.parse(self)

        # Have packages handle user arguments.
        for pkg in self.packages:
            pkg.parse_arguments(self.arguments)

    def argument(self, name):
        val = getattr(self.arguments, name, None)
        if val is None:
            return self._def_args.get(name, None)
        return val

    ##
    ## Search for packages.
    ##
    def configure(self):

        # Only continue if there is a need to reconfigure.
        if not self.needs_configure():
            self.show_configuration(sys.stdout)
            return

        # Check if we need to show configuration and dump.
        self.show_configuration(sys.stdout)

        sys.stdout.write('Configuring...\n')
        sys.stdout.write('  Packages to be configured:\n')
        for pkg in self.packages:
            if pkg.explicit:
                sys.stdout.write('    ' + pkg.name + '\n')

        # Check for downloads.
        sys.stdout.write('  Installing packages...\n')
        for pkg in self.packages:
            pkg.check_download()

        # Repeat the 'search' operation on each package until all return
        # True, indicating there are no new potential locations.
        sys.stdout.write('  Searching for candidates...')
        sys.stdout.flush()
        done = False
        while not done:
            done = True
            for pkg in self.packages:
                if not pkg.search():
                    done = False
        sys.stdout.write(' done.\n')

        # Now we have a list of all locations that passed the footprint
        # testing. Perform thorough checks.
        sys.stdout.write('  Checking candidates...')
        sys.stdout.flush()
        for pkg in self.packages:
            pkg.check()
        sys.stdout.write(' done.\n')

        # We have our list of available packages now. Analyse the entire build
        # graph to determine which packages need to be able to build with each
        # other.
        if self.resolver is not None:
            sys.stdout.write('  Resolving installations...')
            sys.stdout.flush()
            self.resolver(self)
            sys.stdout.write(' done.\n')

        # Clear CRCs and old builders.
        self.crcs = {}
        self.src_crcs = {}
        self.old_bldrs = {}

        # Save configuration results.
        self.save()

        sys.stdout.write('  Success.\n')
        sys.stdout.write('  Configuration details:\n')
        self.write_configuration(sys.stdout, 4)

        # Write directions.
        sys.stdout.write('\n^^^ Scroll up to see the results of configuring. ^^^\n')
        sys.stdout.write('From here you can:\n')
        sys.stdout.write('  Build:              use\n')
        sys.stdout.write('  Show configuration: use -s\n')
        sys.stdout.write('  Reconfigure:        use configure [options]\n')
        sys.stdout.write('  Show help:          use -h\n')

        # Now exit.
        sys.exit()

    ##
    ##
    ##
    def needs_configure(self):

        # Check if the user requested configuration.
        if 'configure' in self.arguments.targets:
            sys.stdout.write('User requested configuration.\n')
            return True

        # If we have not run before then definitely configure.
        if not os.path.exists('.use.db'):
            sys.stdout.write('No prior configuration to use.\n')
            return True

        # Check if anything has changed in the build structure.
        if os.path.exists('.use.db'):
            with open('.use.db', 'r') as inf:
                old_ctx = pickle.load(inf)
        else:
            old_ctx = None
        if self != old_ctx:
            sys.stdout.write('Build structure has changed.\n')
            return True
        self._use_old_ctx(old_ctx)

        # Check each package for reconfiguration requests.
        for pkg in self.packages:
            if pkg.needs_configure():
                sys.stdout.write(pkg.name + ' has changed.\n')
                return True

        return False

    ##
    ## Scan files for rule sources.
    ##
    def find_sources(self):
        for rule in self.rules:
            rule.find_sources(self)

    ##
    ## Expand sources into objects.
    ##
    def expand(self):
        logging.debug('Context: Expanding rules.')
        for rule in self.rules:
            rule.expand(self)

    def scan(self):
        logging.debug('Context: Scanning for dependencies.')
        done = False
        while not done:
            done = True
            for rule in self.rules:
                tmp = rule.scan(self)
                done = tmp if tmp is False else done
        logging.debug('Context: Done scanning for dependencies.')

    ##
    ## After creating all the flows we need to augment them to include
    ## package dependencies.
    ##
    def augment(self):
        logging.debug('Context: Beginning augmentation.')

        # First modify all uses that are contained in a tree.
        for use in self.uses:
            if not isinstance(use, Use):
                continue
            deps = use.package.all_dependencies
            if deps and use.parents:
                logging.debug('Context: Augmenting ' + str(use))
                new_use = use
                parents = list(use.parents)
                for dep in deps:
                    logging.debug('Context: Adding ' + dep.name)

                    # Use existing 'use' if possible.
                    if len(dep.uses):
                        new_use = new_use + dep.uses[0]
                    else:
                        new_use = new_use + Use(dep, None, use.condition)

                for par in parents:
                    if par.left is use:
                        par.left = new_use
                    else:
                        par.right = new_use

        # Now process any uses directly connected to rules. These will not
        # have been picked up in the last run.
        for rule in self.rules:
            use = rule.use
            if not isinstance(use, Use):
                continue
            deps = use.package.dependencies
            if deps:
                logging.debug('Context: Augmenting ' + str(use))
                new_use = use
                for dep in deps:
                    logging.debug('Context: Adding ' + dep.name)

                    # Use existing 'use' if possible.
                    if len(dep.uses):
                        new_use = new_use + dep.uses[0]
                    else:
                        new_use = new_use + Use(dep, None, use.condition)

                rule.use = new_use

        logging.debug('Context: Done augmentation.')

    ##
    ## Decide which targets to build.
    ##
    def find_targets(self):
        logging.debug('Context: Finding targets.')

        # Clear any existing targets.
        self.targets = []

        # Strip out any known dummy targets.
        tgts = set(self.arguments.targets)
        if 'configure' in tgts:
            tgts.remove('configure')

        # Can I find any of the remaining targets in my
        # set of nodes?
        for tgt in tgts:
            node = self.find_node(tgt)
            if node is not None:
                self.targets.append(node)

        # By default we locate any nodes with no children.
        if not self.targets:
            for rule in self.rules:
                for node in rule.product_nodes:
                    if not node.products:
                        self.targets.append(node)
        logging.debug('Context: Found targets: ' + str(self.targets))

        logging.debug('Context: Done finding targets.')

    ##
    ## Build targets.
    ##
    def build(self):
        logging.debug('Context: Building targets.')

        for tgt in self.targets:
            tgt.build(self)

        # Save state.
        self.update_node_crcs()
        self.save()

        logging.debug('Context: Done building targets.')

    def new_arguments(self):
        return Arguments(self)

    def new_options(self, *args, **kwargs):
        return OptionDict(*args, **kwargs)

    def load_package(self, pkg_name, explicit=False):

        # Load the package class.
        pkg_class = load_class(pkg_name)

        # Do we already have this package loaded?
        if pkg_class not in self._pkg_map:

            # Instantiate the package and insert into mapping.
            pkg = pkg_class(self, explicit)
            self._pkg_map[pkg_class] = pkg
            self.packages.append(pkg)

        # Use the existing package.
        else:
            pkg = self._pkg_map[pkg_class]
            pkg.explicit = True if explicit else pkg.explicit

        return pkg

    ##
    ## Create a Use on this context.
    ##
    def new_use(self, pkg_name, opts=None, cond=None, **kwargs):
        pkg = self.load_package(pkg_name, True)
        if opts and kwargs:
            opts = opts + self.new_options(**kwargs)
        elif kwargs:
            opts = self.new_options(**kwargs)
        use = Use(pkg, opts, cond)
        self.uses.append(use)
        return use

    def new_rule(self, src, use, *args, **kwargs):

        # # We must be able to locate the requested use already
        # # in the graph.
        # assert use in self.uses

        # Create a new rule to encapsulate this.
        opts = self.new_options(**kwargs)
        cond = args[0] if len(args) else None
        rule = Rule(src, use, cond, options=opts)

        # # Source can either be a regular expression, a list
        # # of regular expressions or a RuleList. If we have a rule
        # # list we add edges between the rules.
        # if isinstance(src, RuleList):
        #     rule.sub_rules.extend([s for s in src])

        # # Otherwise make sure we have a list and add the rule
        # # as a node.
        # else:
        # rule.source = to_list(src)

        # Store the new rule.
        self.rules.append(rule)
        return rule

    def package(self, name):
        pkg_class = load_class(name)
        return self._pkg_map[pkg_class]

    def exit(self, save=True):
        if save:
            self.update_node_crcs()
            self.save()
        sys.exit(1)

    def node_crc(self, node):
        return self.crcs.get(repr(node), None)

    def node_source_crcs(self, node):
        return self.src_crcs.get(repr(node), None)

    def update_node_crcs(self):
        for n in self._node_map.itervalues():
            self.update_node_crc(n)

        to_del = [k for k, v in self.crcs.iteritems() if v is None]
        for k in to_del:
            del self.crcs[k]

        to_del = [k for k, v in self.src_crcs.iteritems() if v is None]
        for k in to_del:
            del self.src_crcs[k]

    def update_node_crc(self, node):

        # Don't update the node if it has not been seen yet.
        if node.seen:
            self.crcs[repr(node)] = node.current_crc(self)
            self.src_crcs[repr(node)] = node.current_source_crcs(self)

    ##
    ## Store context state to file.
    ##
    def save(self):

        # Parser won't pickle.
        parser = self.parser
        targets = self.arguments.targets
        _arg_map = self._arg_map
        del self.parser
        del self._arg_map
        self.arguments.targets = None

        with open('.use.db', 'w') as out:
            pickle.dump(self, out)

        # Reset.
        self.parser = parser
        self._arg_map = _arg_map
        self.arguments.targets = targets

    ##
    ## Load context from file.
    ##
    def load(self, base_dir):
        pass

    def _use_old_ctx(self, old_ctx):

        # Copy over all the package contents.
        for pkg in self.packages:
            other_pkg = old_ctx._pkg_map[pkg.__class__]
            for ver, other_ver in zip(pkg.versions, other_pkg.versions):
                ver.installations = other_ver.installations
                for inst in ver.installations:

                    # Update version.
                    inst.version = ver

                    # Update features.
                    ftr_names = [f.name for f in inst.features]
                    inst.features = []
                    inst._ftr_map = {}
                    for ftr in ver.features:
                        if ftr.name in ftr_names:
                            inst.features.append(ftr)
                            inst._ftr_map[ftr.name] = ftr

        # Copy over CRCs.
        self.crcs = old_ctx.crcs
        self.src_crcs = old_ctx.src_crcs

        # Build a suite of old builders.
        self.old_bldrs = old_ctx.old_bldrs
        for rule in old_ctx.rules:
            for n in rule.product_nodes:
                if n.builder is not None:
                    self.old_bldrs[repr(n)] = n.builder

        # Run the resolver again.
        if self.resolver is not None:
            self.resolver(self)

    def node(self, node_class, *args, **kwargs):
        n = node_class(*args, **kwargs)
        rep = repr(n)
        return self._node_map.setdefault(rep, n)

    def file(self, *args, **kwargs):
        return self.node(File, *args, **kwargs)

    def write_configuration(self, strm, indent=0):
        for pkg in self.packages:
            insts = list(pkg.iter_installations())
            found = len(insts) > 0
            if found and pkg.explicit:
                strm.write(indent*' ' + '%s\n'%pkg.name)

                # Just use the first installation I can find.
                indent += 2
                if pkg.sub_packages:
                    strm.write(indent*' ' + 'selected package: %s\n'%insts[0].version.package.name)
                if insts[0].location is not None:
                    strm.write(indent*' ' + '%s\n'%insts[0].location.text(indent))
                if insts[0].features:
                    strm.write(indent*' ' + 'features: %s\n'%[f.name for f in insts[0].features])
                indent -= 2

    def show_configuration(self, strm, indent=0):

        # If the user requested to see current configuration then
        # do so now, unless they wished to clean.
        if self.arguments.show_config:
            strm.write(indent*' ' + 'Current arguments:\n')
            indent += 2
            for k, v in self.arguments.__dict__.iteritems():
                if k in ['show_config']:
                    continue
                if v is not None and (not isinstance(v, list) or len(v) > 0):
                    arg = self._arg_map.get(k)
                    strm.write(indent*' ' + '{} {}\n'.format(arg.option_strings[0], v))
            indent -= 2
            strm.write('Current configuration:\n')
            self.write_configuration(strm, indent + 2)
            self.exit(False)

    def norm_path(self, path):
        if os.path.isabs(path):
            cwd = os.getcwd()
            if cwd[-1] != '/':
                cwd += '/'
            com = os.path.commonprefix([cwd, path])
            if com == cwd:
                return os.path.relpath(path, cwd)
        return path

    def find_node(self, path):
        path = self.norm_path(path)
        return self._node_map.get(path, None)
