# cz.pyx
from libc.stdio cimport fopen, fread, fwrite, fclose
from libc.stdlib cimport malloc, free
from libc.string cimport memset

def read_binary_file(str filename, int size):
    cdef FILE *fp
    cdef char *buffer
    cdef int bytes_read
    cdef list result

    fp = fopen(filename.encode('utf-8'), "rb")
    if not fp:
        raise IOError("Failed to open file")

    buffer = <char *>malloc(size)
    if not buffer:
        fclose(fp)
        raise MemoryError("Failed to allocate buffer")

    memset(buffer, 0, size)

    bytes_read = fread(buffer, 1, size, fp)
    result = [buffer[i] for i in range(bytes_read)]

    free(buffer)
    fclose(fp)
    return result

def write_binary_file(str filename, bytes data):
    cdef FILE *fp
    cdef int size = len(data)

    fp = fopen(filename.encode('utf-8'), "wb")
    if not fp:
        raise IOError("Failed to open file for writing")

    fwrite(data, 1, size, fp)
    fclose(fp)
