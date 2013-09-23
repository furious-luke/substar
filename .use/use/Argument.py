import argparse

def boolean(string):
    string = string.lower()
    if string in ['0', 'f', 'false', 'no', 'off']:
        return False
    elif string in ['1', 't', 'true', 'yes', 'on']:
        return True
    else:
        raise ValueError()

class ConfigureAction(argparse.Action):
    
    def __init__(self, option_strings, dest, required=False, help=None, metavar=None, default=None):
        strings = []
        self.positive_strings = set()
        self.negative_strings = set()
        for string in option_strings:
            assert string.startswith('--enable') or string.startswith('--with')
            strings.append(string)
            self.positive_strings.add(string)
            neg_string = string.replace('--enable', '--disable')
            neg_string = neg_string.replace('--with', '--without')
            strings.append(neg_string)
            self.negative_strings.add(neg_string)
        super(ConfigureAction, self).__init__(option_strings=strings, dest=dest, nargs='?', const=None, default=default, type=boolean,
                                              choices=None, required=required, help=help, metavar=metavar)

    def __call__(self, parser, namespace, value, option_string=None):
        if value is None:
            value = option_string in self.positive_strings
        elif option_string in self.negative_strings:
            value = not value
        setattr(namespace, self.dest, value)

class Arguments(object):

    def __init__(self, ctx):
        self._ctx = ctx

    def __call__(self, *args, **kwargs):

        # Replace any default with None and store the default
        # for later.
        default = kwargs.pop('default', None)
        kwargs['default'] = None
        if default is None:
            if kwargs.get('action', None) == 'store_true':
                default = False
            elif kwargs.get('action', None) == 'store_false':
                default = True

        # If the action is 'boolean' then use an augmented boolean type.
        if kwargs.get('action', None) == 'boolean':
            kwargs['action'] = ConfigureAction

        # Create the argument to get hold of destination name.
        act = self._ctx.parser.add_argument(*args, **kwargs)

        # Set the default.
        self._ctx._def_args[act.dest] = default
        self._ctx._arg_map[act.dest] = act

        new_arg = Argument(act.dest, self._ctx)
        setattr(self, new_arg.name, new_arg)
        return self

class Argument(object):

    def __init__(self, name, ctx, use=None, kwargs={}):
        self.name = name
        self.context = ctx
        self.use = use
        self._kwargs = kwargs

    def __eq__(self, op):
        return ArgumentCheck('eq', self, op)

    def __deepcopy__(self, memo):
        return Argument(self.name, self.context, self.use)

    def __add__(self, op):
        return ArgumentCheck('add', self, op)

    def __str__(self):
        return self.value()

    def contains(self, op):
        return ArgumentCheck('in', self, op)

    def value(self):
        if self.use is None:
            return self.context.argument(self.name)
        else:
            attr = getattr(self.use, self.name)
            if callable(attr):
                return attr(**self._kwargs)
            else:
                return attr

    def compare(self, op):
        return self.name == op.name and self.use == op.use

class ArgumentCheck(object):

    def __init__(self, op, left, right=None):
        self.op = op
        self.left = left
        self.right = right

    def __nonzero__(self):
        left_val = self.left.value() if isinstance(self.left, Argument) else self.left
        if self.op == 'eq':
            right_val = self.right.value() if isinstance(self.right, Argument) else self.right
            return left_val == right_val
        elif self.op == 'in':
            right_val = self.right.value() if isinstance(self.right, Argument) else self.right
            return right_val in left_val
        assert 0

    def __eq__(self, op):
        return self.compare(op)

    def __ne__(self, op):
        return not self.compare(op)

    def __str__(self):
        assert self.op in ['add']
        return str(self.value())

    def value(self):
        left_val = self.left.value() if isinstance(self.left, Argument) else self.left
        right_val = self.right.value() if isinstance(self.right, Argument) else self.right
        if self.op == 'add':
            return left_val + right_val
        assert 0

    def compare(self, op):
        if self.op != op.op:
            return False
        if type(self.left) != type(op.left):
            return False
        if isinstance(self.left, (Argument, ArgumentCheck)):
            if not self.left.compare(op.left):
                return False
        elif self.left != op.left:
            return False
        if type(self.right) != type(op.right):
            return False
        if isinstance(self.right, (Argument, ArgumentCheck)):
            if not self.right.compare(op.right):
                return False
        else:
            return self.right == op.right

    def update_context(self, ctx):
        if isinstance(self.left, Argument):
            self.left.context = ctx
        if isinstance(self.right, Argument):
            self.right.context = ctx
