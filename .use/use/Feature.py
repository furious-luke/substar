import logging
from Package import Version
from Use import Use
from .apply import apply

class Feature(Version):

    def __init__(self, version, *args, **kwargs):
        super(Feature, self).__init__(version.package, *args, **kwargs)

        # Swap version name to feature name.
        del self.version
        if not hasattr(self, 'name'):
            name = str(self.__class__)
            self.name = name[name.rfind('.') + 1:name.rfind('\'')].lower()

        # Add the version as an attribute.
        self.version = version

        # No need for patterns.
        del self.patterns

    ##
    ## Check for existence of this feature.
    ##
    def check(self, inst):
        logging.debug('Feature: Checking ' + self.name)
        if inst.location is not None:
            res = self.footprint(inst.location)
        else:
            res = (True,)
        logging.debug('Feature: Done checking ' + self.name)
        return res[0]

    def expand(self, nodes, use_options={}, rule_options={}):
        return self.version.expand(nodes, self, use_options, rule_options)

    def apply(self, prods, use_options={}, rule_options={}):
        apply(self, prods, use_options, rule_options)

    def actions(self, node, target, opts):
        return []

class FeatureUse(Use):

    def __init__(self, name, use, options=None, cond=None):
        self.feature_name = name
        self.use = use
        super(FeatureUse, self).__init__(use.package, options, cond)
