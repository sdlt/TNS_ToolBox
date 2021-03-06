import os
from cffi import FFI
import pathlib

os.system('cp ../libTNS.so ./')
ffi = FFI()
this_dir = pathlib.Path().absolute()

h_file_name = this_dir / "pyTNS.h"
with open(h_file_name) as h_file:
    ffi.cdef(h_file.read())

ffi.set_source(
    "pyTNS",
    '#include "pyTNS.h"',
    libraries=["TNS"],
    library_dirs=[this_dir.as_posix()],
    extra_link_args=["-Wl,-rpath,."])

if __name__ == "__main__":
    ffi.compile(verbose=True)
