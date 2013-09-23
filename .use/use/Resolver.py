import sys
from .Node import Node
from .Use import UseGroup
import logging

##
## Resolve which packages will be used. When searching
## for packages there will likely be several options for
## many. The Resolver's job is to determine which of the
## packages may work in tandem and update all Use objects
## to reflect this selection. This will be a compilcated
## operation, generally, as it will depend on what options
## each use has requested of the same package and also
## what dependencies each package is built against.
##
class Resolver(Node):

    def __init__(self):
        super(Resolver, self).__init__()

    def __eq__(self, op):
        return self.__class__ == op.__class__

    def __ne__(self, op):
        return not self.__eq__(op)

    def __repr__(self):
        return 'resolver'

    ##
    ##
    ##
    def __call__(self, ctx):
        from .Feature import FeatureUse
        logging.debug('Resolver: Resolving package installations.')

        # Build a set of use roots.
        roots = set()
        for rule in ctx.rules:
            self.find_roots(rule.use, roots)

        # Traverse each tree and check
        for root in roots:
            logging.debug('Resolver: Starting at root: ' + str(root))
            if not self.walk(root):
                missing = self.find_missing(root)
                sys.stdout.write('\n    Failed to resolve packages.')
                if len(missing) > 1:
                    sys.stdout.write(' Need at least one of the following:')
                for use in missing:
                    sys.stdout.write('\n      %s\n'%use.package.name)
                sys.exit(1)

        # # Just use the first installation of every package.
        # for use in ctx.uses:
        #     self.check_use(use)

        # # Just use the first installation of every package.
        # for use in ctx.uses:

        #     # Feature uses get handled a little differently.
        #     if isinstance(use, FeatureUse):

        #         # If there is no selected installation then fail.
        #         if use.use.selected is None:
        #             sys.stdout.write('\n    Failed to resolve "' + use.ftr + '".')
        #             sys.exit(1)

        #         # Take the selected installation and see if we can find this
        #         # feature in its list.
        #         use.selected = None
        #         for ftr in use.use.selected.features:
        #             if ftr.name == use.feature_name:
        #                 use.selected = ftr
        #                 break
        #         if use.selected is None:
        #             sys.stdout.write('\n    Failed to resolve "' + use.ftr + '".')
        #             sys.exit(1)

        logging.debug('Resolver: Done resolving package installations.')

    def check_use(self, use):
        from .Feature import FeatureUse
        logging.debug('Resolver: Checking Use.')

        # Don't process feature uses.
        if not isinstance(use, FeatureUse):
            use.selected = self.check_package(use.package)
            use._found = use.selected is not None

    def check_package(self, pkg, fail=True):
        logging.debug('Resolver: Checking package.')

        # Get the first installation.
        sel = None
        for inst in pkg.iter_installations():
            sel = inst
            logging.debug('Resolver: Found installation: ' + str(inst))
            break

        # # Handle no results.
        # if sel is None:
        #     if fail:
        #         sys.stdout.write('\n    Failed to resolve "' + pkg.name + '".\n')
        #         sys.exit(1);

        return sel

    def check_feature(self, use):
        if use.use.enabled:
            for ftr in use.use.selected.features:
                if ftr.name == use.feature_name:
                    use.selected = ftr
                    break

    def find_roots(self, use, roots):
        if not use.parents:
            roots.add(use)
        for par in use.parents:
            self.find_roots(par, roots)

    def walk(self, use):
        from .Feature import FeatureUse
        logging.debug('Resolver: Walking: ' + str(use))
        if isinstance(use, UseGroup):
            logging.debug('Resolver: Have UseGroup.')
            have_left = self.walk(use.left) if use.left is not None else False
            have_right = self.walk(use.right) if use.right is not None else False
            if not have_left:
                if use.op == 'or':
                    use._found = have_right
                    return have_right
                else:
                    use._found = False
                    return False
            if use.op != 'or':
                use._found = have_right
                return have_right
            else:
                use._found = True
                return True
        elif isinstance(use, FeatureUse):
            logging.debug('Resolver: Have FeatureUse.')
            self.check_use(use.use)
            self.check_feature(use)
            return use.selected is not None
        else:
            logging.debug('Resolver: Have Use.')
            self.check_use(use)
            return use.selected is not None

    def find_missing(self, use):
        from .Feature import FeatureUse
        missing = []
        if not isinstance(use, FeatureUse):
            if not use._found:
                if isinstance(use, UseGroup):
                    if use.op == 'or':
                        missing.extend(self.find_missing(use.left))
                        missing.extend(self.find_missing(use.right))
                    else:
                        missing.extend(self.find_missing(use.left))
                        if not missing:
                            missing.extend(self.find_missing(use.right))
                else:
                    missing.append(use)
        else:
            if use.selected is None:
                missing.append(use)
        return missing
