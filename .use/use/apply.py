import os
from Platform import platform

def apply(obj, prods, use_options={}, rule_options={}):
    for prod in prods:
        bldr = prod[1]
        opts = bldr.options
        append_headers(obj, opts)
        if 'compile' not in opts:
            append_libraries(obj, opts)
            append_rpaths(obj, opts)

def append_headers(obj, opts):
    hdr_dirs = getattr(obj, 'header_dirs', [])
    if hdr_dirs:
        dst = opts.setdefault('header_dirs', [])
        for d in hdr_dirs:
            if d not in dst and d not in platform.system_header_dirs:
                dst.append(d)

def append_libraries(obj, opts):
    lib_dirs = getattr(obj, 'library_dirs', [])
    libs = getattr(obj, 'libraries', [])
    if lib_dirs:
        dst = opts.setdefault('library_dirs', [])
        for d in lib_dirs:
            if d not in dst and d not in platform.system_library_dirs:
                dst.append(d)
        opts['library_dirs'] = platform.order_library_dirs(dst)
    if libs:
        dst = opts.setdefault('libraries', [])
        for h in libs:
            if h not in dst:
                dst.append(h)

def append_rpaths(obj, opts):
    lib_dirs = getattr(obj, 'library_dirs', [])
    if lib_dirs:
        dst = opts.setdefault('rpath_dirs', [])
        for d in lib_dirs:
            if d not in dst and d not in platform.system_library_dirs:
                dst.append(d)
        opts['rpath_dirs'] = platform.order_library_dirs(dst)
