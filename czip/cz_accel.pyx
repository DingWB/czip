# cython: language_level=3
import struct
import zlib

# constants used by cz.py
_block_magic = b"MB"

def compress_block(block, level=6):
    """Compress a block (deflate, no headers) and return compressed bytes.

    This mirrors the behaviour in Writer._write_block.
    """
    c = zlib.compressobj(level, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0)
    compressed = c.compress(block) + c.flush()
    return compressed


def load_bcz_block(handle, decompress=False):
    """Read a BCZ block from `handle`.

    If `decompress` is True, return (block_size, data_bytes). Otherwise
    return (block_size, data_len).
    """
    magic = handle.read(2)
    if not magic or magic != _block_magic:
        raise StopIteration
    block_size = struct.unpack("<H", handle.read(2))[0]
    if decompress:
        deflate_size = block_size - 6
        d = zlib.decompressobj(-15)
        data = d.decompress(handle.read(deflate_size)) + d.flush()
        data_len = struct.unpack("<H", handle.read(2))[0]
        return block_size, data
    else:
        # seek forward deflate data and read trailing uncompressed length
        handle.seek(block_size - 6, 1)
        data_len = struct.unpack("<H", handle.read(2))[0]
        return block_size, data_len


def unpack_records(data, fmt):
    """Unpack `data` (bytes) into a list of tuples according to `fmt`.

    `fmt` should be the struct format string without endianness (e.g. 'IIf').
    This uses little-endian ('<') to match the rest of the code.
    """
    if not data:
        return []
    unit = struct.calcsize(fmt)
    n = len(data) // unit
    res = []
    off = 0
    unpack_from = struct.unpack_from
    full_fmt = f"<{fmt}"
    for i in range(n):
        res.append(unpack_from(full_fmt, data, off))
        off += unit
    return res


cdef unsigned long _BLOCK_MAX_LEN = 65535


def c_read(handle, block_raw_length, buffer, within_block_offset, size):
    """Read `size` bytes from BGZF-like stream using existing handle and
    `load_bcz_block` for subsequent blocks. Returns a tuple:
    (data_bytes, new_block_raw_length, new_buffer, new_within_block_offset)
    """
    cdef Py_ssize_t need = size
    cdef bytearray out = bytearray()
    cdef bytes buf = buffer if buffer is not None else b""
    cdef Py_ssize_t buflen = len(buf)
    cdef Py_ssize_t within = within_block_offset
    cdef Py_ssize_t take
    cdef object block_size_data
    while need and block_raw_length:
        buflen = len(buf)
        if within + need <= buflen:
            out.extend(buf[within: within + need])
            within += need
            need = 0
            break
        else:
            # take rest of buffer
            if within < buflen:
                out.extend(buf[within:])
                need -= (buflen - within)
            # load next block
            try:
                block_size_data = load_bcz_block(handle, True)
            except StopIteration:
                # EOF
                buf = b""
                block_raw_length = 0
                within = 0
                break
            block_raw_length, buf = block_size_data
            within = 0

    return bytes(out), block_raw_length, buf, within


def c_readline(handle, block_raw_length, buffer, within_block_offset, newline=b"\n"):
    """Read a single line (including newline) from BGZF-like stream.
    Returns (line_bytes, new_block_raw_length, new_buffer, new_within_block_offset)
    """
    cdef bytearray out = bytearray()
    cdef bytes buf = buffer if buffer is not None else b""
    cdef Py_ssize_t within = within_block_offset
    cdef Py_ssize_t i
    cdef object block_size_data
    cdef Py_ssize_t buflen
    while block_raw_length:
        buflen = len(buf)
        i = buf.find(newline, within)
        if i == -1:
            # append rest and load next block
            if within < buflen:
                out.extend(buf[within:])
            try:
                block_size_data = load_bcz_block(handle, True)
            except StopIteration:
                # EOF
                buf = b""
                block_raw_length = 0
                within = 0
                break
            block_raw_length, buf = block_size_data
            within = 0
        elif i + 1 == buflen:
            # newline at end of block
            out.extend(buf[within:])
            try:
                block_size_data = load_bcz_block(handle, True)
            except StopIteration:
                buf = b""
                block_raw_length = 0
                within = 0
                break
            block_raw_length, buf = block_size_data
            within = 0
            break
        else:
            out.extend(buf[within: i + 1])
            within = i + 1
            break

    return bytes(out), block_raw_length, buf, within


def c_pos2id(handle, block_virtual_offsets, fmts, unit_size, positions, col_to_query=0):
    """Accelerated implementation of pos2id.

    Parameters:
    - handle: open binary file handle positioned anywhere (we'll seek as needed)
    - block_virtual_offsets: iterable of virtual offsets (ints)
    - fmts: format string (like 'IIf')
    - unit_size: record size in bytes
    - positions: iterable of (start, end) pairs
    - col_to_query: index of column to compare

    Returns list where each element is None or [id_start, id_end]
    """
    cdef list results = []
    cdef list block_offsets = list(block_virtual_offsets)
    cdef Py_ssize_t nblocks = len(block_offsets)
    cdef list block_1st_starts = []
    cdef Py_ssize_t idx
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    cdef object block_size_data
    cdef bytes data
    cdef Py_ssize_t off
    cdef object rec
    cdef object val
    cdef Py_ssize_t start_block_index = 0
    cdef Py_ssize_t primary_id
    cdef Py_ssize_t id_start, id_end
    cdef Py_ssize_t i

    # read first record of each block to get starting positions
    for vo in block_offsets:
        block_start = vo >> 16
        within = vo & 0xFFFF
        handle.seek(block_start)
        try:
            block_size_data = load_bcz_block(handle, True)
        except StopIteration:
            block_1st_starts.append(-1)
            continue
        _, data = block_size_data
        if within + unit_size <= len(data):
            rec = struct.unpack(f"<{fmts}", data[within: within + unit_size])
            block_1st_starts.append(rec[col_to_query])
        else:
            # malformed: no record
            block_1st_starts.append(-1)

    # iterate positions
    for start, end in positions:
        # advance start_block_index until next block start > start
        while start_block_index < nblocks - 1 and block_1st_starts[start_block_index + 1] <= start:
            start_block_index += 1
        vo = block_offsets[start_block_index]
        block_start = vo >> 16
        within = vo & 0xFFFF
        # seek to block start and load
        handle.seek(block_start)
        try:
            block_size_data = load_bcz_block(handle, True)
        except StopIteration:
            results.append(None)
            continue
        _, data = block_size_data
        off = within
        # move forward until record[col] >= start
        try:
            rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
        except Exception:
            results.append(None)
            continue
        while rec[col_to_query] < start:
            off += unit_size
            primary_id = (( _BLOCK_MAX_LEN * start_block_index) + off) // unit_size
            if off + unit_size > len(data):
                # load next block
                try:
                    block_size_data = load_bcz_block(handle, True)
                    _, data = block_size_data
                    start_block_index += 1
                    off = 0
                except StopIteration:
                    break
            try:
                rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
            except Exception:
                break
        if rec[col_to_query] < start:
            results.append(None)
            continue
        # compute primary id for start
        if off == 0:
            primary_id = int((_BLOCK_MAX_LEN * start_block_index + off) / unit_size)
        else:
            primary_id = int((_BLOCK_MAX_LEN * start_block_index + off) / unit_size)
        id_start = primary_id
        # advance until record[col] >= end
        while True:
            try:
                off += unit_size
                primary_id += 1
                if off + unit_size > len(data):
                    # load next block
                    try:
                        block_size_data = load_bcz_block(handle, True)
                        _, data = block_size_data
                        start_block_index += 1
                        off = 0
                    except StopIteration:
                        break
                rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
                if rec[col_to_query] >= end:
                    break
            except Exception:
                break
        id_end = primary_id
        results.append([id_start, id_end])

    return results


def c_adjust_virtual_offsets(delta_offset, offsets):
    """Adjust a list of virtual offsets by delta_offset added to the block start."""
    cdef list out = []
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    for vo in offsets:
        block_start = vo >> 16
        within = vo & 0xFFFF
        new_block = block_start + delta_offset
        new_vo = (new_block << 16) | within
        out.append(new_vo)
    return out


def c_pack_records(rows, fmt):
    """Pack a sequence of rows (iterable of sequences) into bytes using struct.pack.

    This loops in C to reduce Python overhead.
    """
    cdef bytearray out = bytearray()
    cdef object row
    cdef object p
    cdef str fullfmt = f"<{fmt}"
    for row in rows:
        # convert to tuple/list to ensure indexability
        try:
            p = struct.pack(fullfmt, *row)
        except Exception:
            # fall back to safer per-element conversion
            p = struct.pack(fullfmt, *tuple(row))
        out.extend(p)
    return bytes(out)


def c_pack_records_fast(rows, fmt):
    """Pack rows into a single bytes object using a pre-allocated buffer

    This reduces per-row allocation by using struct.pack_into on a
    pre-sized bytearray.
    """
    cdef object rows_list = list(rows)
    if not rows_list:
        return b""
    fullfmt = f"<{fmt}"
    unit = struct.calcsize(fullfmt)
    n = len(rows_list)
    out = bytearray(unit * n)
    cdef Py_ssize_t i
    for i in range(n):
        try:
            struct.pack_into(fullfmt, out, i * unit, *rows_list[i])
        except Exception:
            struct.pack_into(fullfmt, out, i * unit, *tuple(rows_list[i]))
    return bytes(out)


def c_block_first_values(handle, block_virtual_offsets, fmts, unit_size, s):
    """Return list of first record value at column `s` for each block.

    This seeks to each block, loads it (decompressed), and reads the value.
    """
    cdef list out = []
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    cdef object block_size_data
    cdef bytes data
    cdef object rec
    for vo in block_virtual_offsets:
        block_start = vo >> 16
        within = vo & 0xFFFF
        handle.seek(block_start)
        try:
            block_size_data = load_bcz_block(handle, True)
        except StopIteration:
            out.append(None)
            continue
        _, data = block_size_data
        if within + unit_size <= len(data):
            rec = struct.unpack(f"<{fmts}", data[within: within + unit_size])
            out.append(rec[s])
        else:
            out.append(None)
    return out


def c_read_1record(handle, block_raw_length, buffer, within, fmts, unit_size):
    """Read a single record (unpacked tuple) from the stream, updating buffer state.

    Returns (record_tuple, new_block_raw_length, new_buffer, new_within)
    """
    cdef bytes buf = buffer if buffer is not None else b""
    cdef Py_ssize_t buflen = len(buf)
    if within + unit_size <= buflen:
        rec = struct.unpack(f"<{fmts}", buf[within: within + unit_size])
        within += unit_size
        return rec, block_raw_length, buf, within
    # need to load next block or remaining bytes
    out = bytearray()
    # take remaining
    if within < buflen:
        out.extend(buf[within:])
    try:
        block_size_data = load_bcz_block(handle, True)
    except StopIteration:
        # EOF
        return None, 0, b"", 0
    block_raw_length, buf = block_size_data
    # fill remaining needed bytes
    needed = unit_size - len(out)
    if needed <= len(buf):
        out.extend(buf[:needed])
        within = needed
        rec = struct.unpack(f"<{fmts}", bytes(out))
        return rec, block_raw_length, buf, within
    else:
        # malformed, not enough
        return None, block_raw_length, buf, 0


def c_seek_and_read_1record(handle, virtual_offset, fmts, unit_size):
    """Seek to virtual offset and read one record (unpacked tuple)."""
    start = virtual_offset >> 16
    within = virtual_offset & 0xFFFF
    handle.seek(start)
    try:
        block_size_data = load_bcz_block(handle, True)
    except StopIteration:
        return None
    _, buf = block_size_data
    if within + unit_size <= len(buf):
        rec = struct.unpack(f"<{fmts}", buf[within: within + unit_size])
        return rec
    # need to assemble across block boundary
    out = bytearray()
    if within < len(buf):
        out.extend(buf[within:])
    try:
        block_size_data = load_bcz_block(handle, True)
    except StopIteration:
        return None
    _, buf2 = block_size_data
    needed = unit_size - len(out)
    if needed <= len(buf2):
        out.extend(buf2[:needed])
        return struct.unpack(f"<{fmts}", bytes(out))
    return None


def c_query_regions(handle, block_virtual_offsets, fmts, unit_size, regions, s, e, dim):
    """Accelerated query for a single dim.

    Returns a list of results where the first matching item for each region
    is ("primary_id_&_dim:", primary_id, dim) followed by (dim, record_tuple)
    for each record whose column[s]..[e] falls into the region.
    """
    cdef list results = []
    cdef list block_offsets = list(block_virtual_offsets)
    cdef Py_ssize_t nblocks = len(block_offsets)
    cdef list block_1st_starts = c_block_first_values(handle, block_offsets, fmts, unit_size, s)
    cdef Py_ssize_t start_block_index = 0
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    cdef object block_size_data
    cdef bytes data
    cdef Py_ssize_t off
    cdef object rec
    cdef Py_ssize_t primary_id
    for start, end in regions:
        # find start block index
        while start_block_index < nblocks - 1 and block_1st_starts[start_block_index + 1] <= start:
            start_block_index += 1
        vo = block_offsets[start_block_index]
        block_start = vo >> 16
        within = vo & 0xFFFF
        handle.seek(block_start)
        try:
            block_size_data = load_bcz_block(handle, True)
        except StopIteration:
            continue
        _, data = block_size_data
        off = within
        # advance to first record >= start
        try:
            rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
        except Exception:
            continue
        while rec[s] < start:
            off += unit_size
            if off + unit_size > len(data):
                # load next block
                try:
                    block_size_data = load_bcz_block(handle, True)
                    _, data = block_size_data
                    start_block_index += 1
                    off = 0
                except StopIteration:
                    break
            try:
                rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
            except Exception:
                break
        if rec[s] < start:
            # nothing found
            continue
        # compute primary id
        primary_id = int((_BLOCK_MAX_LEN * start_block_index + off) / unit_size)
        results.append(("primary_id_&_dim:", primary_id, dim))
        # yield records while record[e] <= end
        while True:
            if rec[e] <= end:
                results.append((dim, rec))
                # move to next record
                off += unit_size
                primary_id += 1
                if off + unit_size > len(data):
                    try:
                        block_size_data = load_bcz_block(handle, True)
                        _, data = block_size_data
                        start_block_index += 1
                        off = 0
                    except StopIteration:
                        break
                try:
                    rec = struct.unpack(f"<{fmts}", data[off: off + unit_size])
                except Exception:
                    break
            else:
                break
    return results


def c_write_chunk_tail(handle, chunk_data_len, block_1st_record_virtual_offsets, chunk_dims):
    """Write chunk tail to `handle`:
    - chunk_data_len (Q)
    - n_blocks (Q)
    - list of block start offsets (n_blocks * Q)
    - list of dims: for each dim, write len as B then bytes(dim)

    This mirrors Writer._write_chunk_tail but performs the packing in Cython.
    """
    # write chunk_data_len and n_blocks
    handle.write(struct.pack("<Q", chunk_data_len))
    n_blocks = len(block_1st_record_virtual_offsets)
    handle.write(struct.pack("<Q", n_blocks))
    # write block offsets
    for vo in block_1st_record_virtual_offsets:
        handle.write(struct.pack("<Q", vo))
    # write dims
    for dim in chunk_dims:
        dim_b = bytes(dim, 'utf-8')
        dim_len = len(dim_b)
        handle.write(struct.pack("<B", dim_len))
        handle.write(struct.pack(f"<{dim_len}s", dim_b))
    return None


# -- Use zlib C API for raw deflate/inflate --
cdef extern from "zlib.h":
    ctypedef unsigned int uInt
    ctypedef unsigned long uLong
    ctypedef void* voidpf

    cdef struct z_stream_s:
        unsigned char *next_in
        uInt avail_in
        uLong total_in
        unsigned char *next_out
        uInt avail_out
        uLong total_out
        char *msg
        voidpf state
        voidpf zalloc
        voidpf zfree
        voidpf opaque

    int deflateInit2_(z_stream_s *strm, int level, int method, int windowBits, int memLevel, int strategy, const char *version, int stream_size)
    int deflate(z_stream_s *strm, int flush)
    int deflateEnd(z_stream_s *strm)
    int inflateInit2_(z_stream_s *strm, int windowBits, const char *version, int stream_size)
    int inflate(z_stream_s *strm, int flush)
    int inflateEnd(z_stream_s *strm)
    const char * zlibVersion()

    int compress2(unsigned char *dest, uLong *destLen, const unsigned char *source, uLong sourceLen, int level)


from cpython.bytearray cimport PyByteArray_FromStringAndSize, PyByteArray_AsString
from cpython.bytes cimport PyBytes_FromStringAndSize
from libc.string cimport memset


def c_inflate_bytes(data):
    """Decompress raw deflate bytes using zlib C inflate with -MAX_WBITS.

    This will attempt to inflate the entire input in a single-pass using a
    conservative output buffer that grows if needed.
    """
    cdef z_stream_s strm
    cdef int ret
    cdef const char *in_ptr = data
    cdef Py_ssize_t in_len = len(data)
    cdef Py_ssize_t out_size = max(in_len * 4 + 1024, 65536)
    cdef object out_ba = PyByteArray_FromStringAndSize(NULL, out_size)
    cdef unsigned char *out_ptr = <unsigned char *> PyByteArray_AsString(out_ba)

    memset(&strm, 0, sizeof(strm))
    strm.next_in = <unsigned char *> in_ptr
    strm.avail_in = <uInt> in_len
    strm.next_out = out_ptr
    strm.avail_out = <uInt> out_size

    ret = inflateInit2_(&strm, -15, zlibVersion(), sizeof(strm))
    if ret != 0:
        raise RuntimeError('inflateInit2_ failed')

    ret = inflate(&strm, 4)  # Z_FINISH == 4
    if ret not in (0, 1):  # Z_OK(0) or Z_STREAM_END(1)
        inflateEnd(&strm)
        raise RuntimeError(f'inflate failed {ret}')

    out_len = out_size - strm.avail_out
    # shrink to actual size
    res = PyBytes_FromStringAndSize(<char *> out_ptr, out_len)
    inflateEnd(&strm)
    return res


def c_deflate_raw(block, level=6):
    """Compress bytes into raw deflate using zlib C deflateInit2_/deflate.

    This attempts a single-pass deflate with a preallocated output buffer.
    """
    cdef z_stream_s strm
    cdef int ret
    cdef const char *in_ptr = block
    cdef Py_ssize_t in_len = len(block)
    cdef Py_ssize_t out_size = in_len + (in_len >> 2) + 1024
    cdef object out_ba = PyByteArray_FromStringAndSize(NULL, out_size)
    cdef unsigned char *out_ptr = <unsigned char *> PyByteArray_AsString(out_ba)

    memset(&strm, 0, sizeof(strm))
    strm.next_in = <unsigned char *> in_ptr
    strm.avail_in = <uInt> in_len
    strm.next_out = out_ptr
    strm.avail_out = <uInt> out_size

    ret = deflateInit2_(&strm, level, 8, -15, 8, 0, zlibVersion(), sizeof(strm))
    if ret != 0:
        raise RuntimeError('deflateInit2_ failed')

    ret = deflate(&strm, 4)  # Z_FINISH
    if ret not in (0, 1):
        deflateEnd(&strm)
        raise RuntimeError(f'deflate failed {ret}')

    out_len = out_size - strm.avail_out
    res = PyBytes_FromStringAndSize(<char *> out_ptr, out_len)
    deflateEnd(&strm)
    return res


def c_fetch_chunk(handle, chunk_start_offset_plus10, block_virtual_offsets, fmts, unit_size):
    """Read all blocks for a chunk and unpack into list of tuples according to fmt.

    Parameters:
    - handle: binary file handle
    - chunk_start_offset_plus10: the offset to the first block (chunk_start + 10)
    - block_virtual_offsets: list of virtual offsets for block first records
    - fmts: struct formats string
    - unit_size: size of each record

    Returns list of tuples (records)
    """
    cdef bytes all_data = b""
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    cdef object block_size_data
    for vo in block_virtual_offsets:
        block_start = vo >> 16
        handle.seek(block_start)
        try:
            block_size_data = load_bcz_block(handle, True)
        except StopIteration:
            break
        _, data = block_size_data
        all_data += data
    # return raw packed data (caller can unpack efficiently)
    if not all_data:
        return b""
    return all_data


def c_get_records_by_ids(handle, chunk_block_virtual_offsets, unit_size, IDs):
    """Fetch records (raw bytes) for given primary IDs. IDs is an iterable of ints (1-based).

    Returns list of bytes objects (each is a record of unit_size bytes).
    """
    cdef list results = []
    cdef Py_ssize_t block_index
    cdef unsigned long vo
    cdef unsigned long block_start
    cdef unsigned long within
    cdef object block_size_data
    cdef bytes buf = b""
    cdef Py_ssize_t current_block_start = -1
    cdef Py_ssize_t offset_in_block
    for ID in IDs:
        idx = ((ID - 1) * unit_size) // _BLOCK_MAX_LEN
        within = ((ID - 1) * unit_size) % _BLOCK_MAX_LEN
        vo = chunk_block_virtual_offsets[idx]
        block_start = vo >> 16
        if block_start != current_block_start:
            handle.seek(block_start)
            try:
                block_size_data = load_bcz_block(handle, True)
            except StopIteration:
                results.append(b"")
                current_block_start = block_start
                buf = b""
                continue
            _, buf = block_size_data
            current_block_start = block_start
        # extract record bytes
        if within + unit_size <= len(buf):
            results.append(buf[within: within + unit_size])
        else:
            # malformed: not enough bytes in this block
            results.append(b"")
    return results
