from cffi import FFI
import pathlib

ffi = FFI()
this_dir = pathlib.Path().absolute()
lib_dir = this_dir / "../"

h_file_name = this_dir / "pyTNS.h"
with open(h_file_name) as h_file:
    ffi.cdef(h_file.read())

ffi.set_source(
    "pyTNS",
    # Since you're calling a fully-built library directly, no custom source
    # is necessary. You need to include the .h files, though, because behind
    # the scenes cffi generates a .c file that contains a Python-friendly
    # wrapper around each of the functions.
    '#include "pyTNS.h"',
    # The important thing is to include the pre-built lib in the list of
    # libraries you're linking against:
    libraries=["TNS"],
    library_dirs=[this_dir.as_posix()],
    extra_link_args=["-Wl,-rpath,."],
)

if __name__ == "__main__":
    ffi.compile(verbose=True)
