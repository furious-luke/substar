import copy, os
from collections import OrderedDict
from File import File
from .Argument import ArgumentCheck, Argument
from .conv import to_list

##
##
##
class Option(object):

    def __init__(self, name=None, short_opts=None, long_opts=None, **kwargs):
        self.name = name
        self.short_opts = to_list(short_opts)
        self.long_opts = to_list(long_opts)
        if self.name is None:
            if self.short_opts:
                self.name = self.short_opts[0]
            else:
                self.name = self.long_opts[0]
        self.space = kwargs.get('space', True)
        self.text = kwargs.get('text', None)
        self.abspath = kwargs.get('abspath', False)

    def __eq__(self, op):
        if isinstance(op, Option):
            return self.name == op.name
        else:
            return self.name == op or op in self.short_opts or op in self.long_opts

    def __hash__(self):
        return hash(self.name)

    def __call__(self, value, opts, short=True):
        opt = self._opt(short)
        str_list = []
        if isinstance(value, bool):
            if value:
                if self.text is not None:
                    str_list.append(self.text.format(**opts))
                else:
                    str_list.append(opt)
        elif isinstance(value, (list, tuple)):
            for v in value:
                self._append(str_list, opt, v)
        else:
            self._append(str_list, opt, value)
        return str_list

    def _append(self, str_list, opt, val):
        if self.abspath:
            if not os.path.isabs(val):
                val = os.path.join(os.getcwd(), val)
        if opt:
            if not self.space:
                str_list.append(opt + str(val))
            else:
                str_list.append(opt)
        if self.space:
            str_list.append(str(val))

    def _opt(self, short):
        if short and self.short_opts:
            return self.short_opts[0]
        elif self.long_opts:
            return self.long_opts[0]
        return None

##
##
##
class OptionParser(object):

    def __init__(self, binary=''):
        self.binary = binary
        self._opts = OrderedDict()

    def add(self, opt):
        self._opts[opt.name] = opt

    def __call__(self, *args, **kwargs):

        # Start with the binary name.
        str_list = [self.binary]

        # The user may have given us a string to process
        # instead of an option dictionary.
        options = {}
        if len(args):
            options = self.make_options_dict(args[0])

        # Insert keyword options.
        options.update(kwargs)

        # Process each option.
        for name, opt in self._opts.iteritems():
            val = options.get(name, None)
            if val is not None:
                str_list.extend(opt(val, options))

        return ' '.join(str_list)

    def make_options_dict(self, opt_str):
        if isinstance(opt_str, (OptionDict, OptionJoin)):
            return opt_str.get()
        elif isinstance(opt_str, dict):
            return dict(opt_str)
        elif opt_str is None:
            return {}
        assert 0, 'Not implemented yet.'

class OptionDict(object):

    def __init__(self, cond=True, **kwargs):
        self.condition = cond
        self._opts = kwargs

    def __repr__(self):
        return str(self.get())

    def __eq__(self, op):
        if type(self) != type(op):
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

        return self._opts == op._opts;

    def __ne__(self, op):
        return not self.__eq__(op)

    def __add__(self, op):
        return OptionJoin(self, op)

    def get(self):
        if bool(self.condition):
            return copy.deepcopy(self._opts)
        else:
            return {}

    def parse(self, ctx):
        for k, v in self._opts.iteritems():
            if isinstance(v, (Argument, ArgumentCheck)):
                try:
                    self._opts[k] = str(v)
                except:
                    import pdb
                    pdb.set_trace()
                    print 'hello'

class OptionJoin(object):

    def __init__(self, left, right):
        self.left = left
        self.right = right
        self._merged = None

    def __repr__(self):
        return str(self.get())

    def __eq__(self, op):
        if type(self) != type(op):
            return False
        if self.left != op.left:
            return False
        if self.right != op.right:
            return False
        return True

    def __ne__(self, op):
        return not self.__eq__(op)

    def __add__(self, op):
        return OptionJoin(self, op)

    def get(self):
        if self._merged is None:
            self._merged = {}
            self._update(self.left)
            self._update(self.right)
            merged = self._merged
            self._merged = None
        # return copy.deepcopy(self._merged)
        return merged

    def parse(self, ctx):
        self.left.parse(ctx)
        self.right.parse(ctx)

    def _update(self, op):
        for k, v in op.get().iteritems():
            if k in self._merged:
                my_val = self._merged[k]
                if isinstance(my_val, list):
                    for x in to_list(v):
                        my_val.append(x)
                else:
                    self._merged[k] = v
            else:
                self._merged[k] = copy.deepcopy(v)
