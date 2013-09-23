import logging
from .Validatable import Validatable

class Node(Validatable):

    def __init__(self, *args, **kwargs):
        super(Node, self).__init__()
        self.rule = None
        self.builder = None
        self.products = []
        self.dependencies = []
        self.scanner = None
        self.seen = False
        self._invalid = False
        self._src_crcs = None
        self._done_scan = False

    def __eq__(self, op):
        return repr(self) == repr(op)

    def __ne__(self, op):
        return not self.__eq__(op)

    def __repr__(self):
        return 'Node'

    ##
    ## Called to process this node.
    ##
    def build(self, ctx):
        logging.debug('Node: Building node: ' + str(self))

        # If we've already processed this node don't do
        # so again.
        if self.seen:
            logging.debug('Node: Already seen this node.')
            return self._invalid

        # Call out to our parents first.
        if self.builder:
            self._invalid = self.builder.build_sources(ctx)
            if self._invalid:
                logging.debug('Node: Parents are invalidated.')

        # Build our dependencies.
        invalid = self.build_dependencies(ctx)
        if invalid:
            self._invalid = invalid
            logging.debug('Node: Dependencies are invalidated.')

        # If our parents are invalidated we must rebuild. If not,
        # check if we're invalidated in any other way.
        if not self._invalid:
            self._invalid = self.invalidated(ctx)
            logging.debug('Node: Set invalidated state to %s'%self._invalid)

        # If we are invalid perform an update.
        if self._invalid:
            self.update(ctx)

        # Flag as seen.
        self.seen = True

        # Return our invalidation state.
        logging.debug('Node: Done building node: ' + str(self))
        return self._invalid

    def build_dependencies(self, ctx):
        logging.debug('Node: Building dependencies.')

        invalid = False
        for dep in self.dependencies:
            this = dep.build(ctx)
            invalid = invalid or this

        logging.debug('Node: Done building dependencies.')
        return invalid

    ##
    ## Determine if this node is invalidated.
    ##
    def invalidated(self, ctx):

        # If any of my sources are invalidated then I must be so. However
        # this is already checked in build.

        # Compare CRCs of sources to those I have stored.
        if self.builder:
            old_src_crcs = ctx.node_source_crcs(self)
            if old_src_crcs is None:
                return True
            for src in self.builder.dependent_nodes:
                if src.builder is None:
                    crc = old_src_crcs.get(repr(src), None)
                    if crc is None or crc != src.current_crc(ctx):
                        return True

            # Also check if the builder has changed.
            if hasattr(ctx, 'old_bldrs'):
                old_bldr = ctx.old_bldrs.get(repr(self), None)
                if old_bldr is None or self.builder != old_bldr:
                    return True

        return False

    ##
    ## Do what is needed to validate this node.
    ##
    def update(self, ctx):
        logging.debug('Node: Updating node: ' + str(self))

        if self.builder:
            self.builder.update(ctx)

        logging.debug('Node: Done updating node: ' + str(self))

    def scan(self, ctx, bldr):
        logging.debug('Node: Scanning.')
        if self.scanner is None:
            scanner = bldr.options.get('scanner', None)
            if scanner is not None:
                scanner = scanner(ctx)
        else:
            scanner = self.scanner
        if not self._done_scan:
            if scanner is not None:
                logging.debug('Node: Using scanner: ' + str(scanner.__class__))
                with open(str(self), 'r') as src_file:
                    data = src_file.read()
                new_deps = list(scanner.find_all(self, data, bldr))
                logging.debug('Node: New dependencies: ' + str(new_deps))
                self.dependencies.extend(new_deps)
                self._new_crc = self._crc32(data)
        logging.debug('Node: Done scanning.')

    def update_source_crcs(self, ctx):
        if self.builder:
            self._src_crcs = {}
            for src in self.builder.dependent_nodes:
                if src.builder is None:
                    self._src_crcs[repr(src)] = src.current_crc(ctx)
        else:
            self._src_crcs = None

    def current_source_crcs(self, ctx):
        if self._src_crcs is None:
            self.update_source_crcs(ctx)
        return self._src_crcs

class Always(Node):

    def __init__(self, *args, **kwargs):
        super(Always, self).__init__(*args, **kwargs)
        self.path = None

    def invalidated(self, ctx):
        return True
