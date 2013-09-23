def find(ctr, match, default=None):
    def _find(ctr, match):
        if hasattr(ctr, '__iter__'):
            if match in ctr:
                if hasattr(ctr, 'itervalues'):
                    # Dereference to value.
                    return (ctr[match], True)
                else:
                    # Return actual object.
                    return (ctr[ctr.index(match)], True)
            if hasattr(ctr, '__getitem__') and hasattr(ctr, 'itervalues'):
                # dict-like object.
                for value in ctr.itervalues():
                    result = _find(value, match)
                    if result[1]:
                        return result
            else:
                # Any other iterable object.
                for value in ctr:
                    result = _find(value, match)
                    if result[1]:
                        return result
        return (None, False)

    result = _find(ctr, match)
    if result[1]:
        return result[0]
    else:
        return default
