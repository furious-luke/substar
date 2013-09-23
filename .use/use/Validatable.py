import os, binascii

class Validatable(object):

    def __init__(self):
        self._crc = None
        self._new_crc = None

    def set_valid_crc(self, filename):
        self._crc = self._crc32_file(filename)

    def invalidated_crc(self, filename, ctx):
        if self._crc is None:
            self._crc = ctx.node_crc(self)
        if self._new_crc is None:
            self._new_crc = self._crc32_file(filename)
        return self._crc is None or self._new_crc is None or self._new_crc != self._crc

    def current_crc(self, ctx):
        return self._new_crc if self._new_crc is not None else self._crc

    def _crc32_file(self, filename):
        if os.path.exists(filename):
            return binascii.crc32(open(filename, 'rb').read()) & 0xFFFFFFFF
        else:
            return None

    def _crc32(self, data):
        return binascii.crc32(data) & 0xFFFFFFFF
