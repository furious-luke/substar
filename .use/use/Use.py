import copy
from .Node import Node
from .Argument import ArgumentCheck, Argument
from .utils import load_class, getarg

##
## Represents a package installation with options.
##
class Use(Node):

    def __init__(self, package, options=None, cond=None):
        super(Use, self).__init__()
        self.package = package
        self.condition = cond
        self.options = options
        self.selected = None
        self.parents = []
        self._found = False
        self.package.uses.append(self)

    def __eq__(self, op):
        if self.package != op.package:
            return False

        # TODO: Fix this condition compare.
        if type(self.condition) != type(op.condition):
            return False
        if isinstance(self.condition, ArgumentCheck):
            if not self.condition.compare(op.condition):
                return False
        elif isinstance(op.condition, ArgumentCheck):
            if not op.condition.coimpare(self.condition):
                return False
        elif self.condition != op.condition:
            return False

        if self.options != op.options:
            return False
        return True

    def __ne__(self, op):
        return not self.__eq__(op)

    def __repr__(self):
        text = 'use<' + repr(self.package)
        if self.options:
            text += ', ' + repr(self.options)
        text += '>'
        return text

    def __add__(self, op):
        grp = UseGroup(self, op, 'add')
        self.parents.append(grp)
        op.parents.append(grp)
        return grp

    def __and__(self, op):
        grp = UseGroup(self, op, 'and')
        self.parents.append(grp)
        op.parents.append(grp)
        return grp

    def __or__(self, op):
        grp = UseGroup(self, op, 'or')
        self.parents.append(grp)
        op.parents.append(grp)
        return grp

    @property
    def found(self):
        return self.selected is not None

    @property
    def enabled(self):
        return self.found and (self.condition is None or bool(self.condition))

    @property
    def have(self):
        return Argument('enabled', self.package.ctx, self)

    def _has_feature(self, name):
        if self.enabled:
            for ftr in self.selected.features:
                if ftr.name == name:
                    return True
        return False

    def has_feature(self, name):
        return Argument('_has_feature', self.package.ctx, self, {'name': name})

    def has_packages(self):
        return True

    def package_iter(self):
        yield self.package

    ##
    ## Make a new feature use.
    ##
    def feature(self, name, opts=None, cond=None, **kwargs):
        if opts and kwargs:
            opts = opts + self.package.ctx.new_options(**kwargs)
        elif kwargs:
            opts = self.package.ctx.new_options(**kwargs)
        return self.package.feature_use(name, self, opts, cond)

    def make_options_dict(self, options):
        return self.selected.options().make_options_dict(options)

    def apply(self, prods, options={}):
        if self.condition is not None and not bool(self.condition):
            return None
        self.selected.apply(prods, self.options, options)
        self.apply_features(prods, options)

    def expand(self, nodes, options={}):
        if self.condition is not None and not bool(self.condition):
            return None
        prods = self.selected.expand(nodes, self.options, options)
        if prods is not None:
            self.apply_features(prods, options)
        return prods

    def apply_features(self, prods, rule_options={}):
        for ftr in self.selected.features:
            ftr.apply(prods, self.options, rule_options)

##
##
##
class UseGroup(object):

    def __init__(self, left, right, op):
        self.left = left
        self.right = right
        self.op = op
        self.parents = []
        self._found = False

    def __eq__(self, op):
        if type(self) != type(op):
            return False
        return self.left == op.left and self.right == op.right and self.op == op.op

    def __ne__(self, op):
        return not self.__eq__(op)

    def __add__(self, op):
        grp = UseGroup(self, op, 'add')
        self.parents.append(grp)
        return grp

    def __and__(self, op):
        grp = UseGroup(self, op, 'and')
        self.parents.append(grp)
        return grp

    def __or__(self, op):
        grp = UseGroup(self, op, 'or')
        self.parents.append(grp)
        return grp

    @property
    def found(self):
        if self.op == 'add':
            return self.left.found and self.right.found
        elif self.op in ['or', 'and']:
            return self.left.found or self.right.found
        return False

    @property
    def enabled(self):
        if self.op == 'add':
            return self.left.enabled and self.right.enabled
        elif self.op in ['or', 'and']:
            return self.left.enabled or self.right.enabled
        return False

    def apply(self, prods):
        if self.op == 'add':
            self.left.apply(prods)
            self.right.apply(prods)

        elif self.op == 'or':
            if self.left.enabled:
                self.left.apply(prods)
            else:
                self.right.apply(prods)

        elif self.op == 'and':
            if self.left.enabled:
                self.left.apply(prods)
            if self.right.enabled:
                if self.left.enabled:
                    right_prods = copy.deepcopy(prods)
                    self.right.apply(right_prods)
                    prods.extend(right_prods)
                else:
                    self.right.apply(prods)

    def expand(self, nodes, options):
        if self.op == 'add':
            prods = self.left.expand(nodes, options)
            if prods is None:
                prods = self.right.expand(nodes, options)
                if prods is not None:
                    self.left.apply(prods)
            else:
                self.right.apply(prods)

        elif self.op == 'or':
            if self.left.enabled:
                prods = self.left.expand(nodes, options)
            else:
                prods = self.right.expand(nodes, options)

        elif self.op == 'and':
            left_prods = self.left.expand(nodes, options) if self.left.enabled else []
            right_prods = self.right.expand(nodes, options) if self.right.enabled else []
            prods = (left_prods if left_prods is not None else []) + (right_prods if right_prods is not None else [])

        return prods
