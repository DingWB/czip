#!/usr/bin/env python
import glob
import os, sys
import struct
import zlib
from builtins import open as _open
import numpy as np
import pandas as pd
import gzip
import math

_bcz_magic = b'BMZIP'
_block_magic = b"MB"
_chunk_magic = b"MC"
_BLOCK_MAX_LEN = 65535
_bcz_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BM\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_version = 1.0


# ==========================================================
def str2byte(x):
    return bytes(str(x), 'utf-8')
# ==========================================================
dtype_func = {
    'h': int, 'H': int,
    'i': int, 'I': int, 'b': int, 'B': int,
    'L': int, 'l': int, 'q': int, 'Q': int,
    'f': float, 'd': float,
    's': str2byte, 'c': str2byte
}


# ==========================================================
def get_dtfuncs(formats,tobytes=True):
    D=dtype_func
    if not tobytes:
        D['s']=str
        D['c']=str
    return [D[t[-1]] for t in formats]


# ==========================================================
def make_virtual_offset(block_start_offset, within_block_offset):
    if within_block_offset < 0 or within_block_offset >= _BLOCK_MAX_LEN:
        raise ValueError(
            "Require 0 <= within_block_offset < 2**16, got %i" % within_block_offset
        )
    if block_start_offset < 0 or block_start_offset >= 2**48:
        raise ValueError(
            "Require 0 <= block_start_offset < 2**48, got %i" % block_start_offset
        )
    return (block_start_offset << 16) | within_block_offset


# ==========================================================
def split_virtual_offset(virtual_offset):
    start = virtual_offset >> 16
    return start, virtual_offset ^ (start << 16)


# ==========================================================
def SummaryBczBlocks(handle):
    if isinstance(handle, Reader):
        raise TypeError("Function BczBlocks expects a binary handle")
    data_start = 0
    while True:
        start_offset = handle.tell()
        try:
            block_length, data_len = _load_bcz_block(handle)
        except StopIteration:
            break
        yield start_offset, block_length, data_start, data_len
        data_start += data_len

# ==========================================================
def _load_bcz_block(handle, decompress=False):
    magic = handle.read(2)
    if not magic or magic != _block_magic:  # next chunk or EOF
        raise StopIteration
    block_size = struct.unpack("<H", handle.read(2))[0]
    # block size is the size of compressed data + 4 (header + tail)
    """
	2 bytes is the size of block_size
	10 bytes for the GZIP header.
	2 bytes for the block tail: block_data_len
	1 byte for the empty GZIP trailer.
	"""
    if decompress:
        """
		# there is no gzip header for deflate compressed format
		to (de-)compress deflate format, use wbits = -zlib.MAX_WBITS (-15)
		to (de-)compress zlib format, use wbits = zlib.MAX_WBITS
		to (de-)compress gzip format, use wbits = zlib.MAX_WBITS | 16
		"""
        deflate_size = block_size - 6
        d = zlib.decompressobj(-15)  # -zlib.MAX_WBITS, means no headers
        data = d.decompress(handle.read(deflate_size)) + d.flush()
        data_len = struct.unpack("<H", handle.read(2))[0]
        # refer to: http://samtools.github.io/hts-specs/SAMv1.pdf
        # data_len = len(data) #uncompressed data length
        return block_size, data
    else:
        handle.seek(block_size - 6, 1)
        data_len = struct.unpack("<H", handle.read(2))[0]
        return block_size, data_len
# ==========================================================
def open1(infile):
    if hasattr(infile, 'readline'):
        return infile
    if infile.endswith('.gz'):
        f = gzip.open(infile, 'rb')
    else:
        f = open(infile, 'r')
    return f
# ==========================================================
def _gz_input_parser(infile,formats,sep='\t',usecols=[1,4,5],dim_cols=[0],
                 chunksize=5000):
    f = open1(infile)
    data, i, pre_dims = [], 0, None
    dtfuncs = get_dtfuncs(formats)
    line = f.readline()
    while line:
        line = line.decode('utf-8')
        values = line.rstrip('\n').split(sep)
        dims = [values[i] for i in dim_cols]
        if dims != pre_dims:  # a new dim (chrom), for example, chr1 -> chr2
            if len(data) > 0:  # write rest data of chr1
                yield pd.DataFrame(data,columns=usecols), pre_dims
                data, i = [], 0
            pre_dims = dims
        if i >= chunksize:  # dims are the same, but reach chunksize
            yield pd.DataFrame(data,columns=usecols), pre_dims
            data, i = [], 0
        data.append([func(v) for v, func in zip([values[i] for i in usecols],
                                                dtfuncs)])
        line = f.readline()
        i += 1
    f.close()
    if len(data)>0:
        yield pd.DataFrame(data, columns=usecols), pre_dims

# ==========================================================
def _text_input_parser(infile,formats,sep='\t',usecols=[1,4,5],dim_cols=[0],
                 chunksize=5000):
    """
    Parse text input file (txt, csv, tsv and so on.) into chunks, every chunk has
    the same dim_cols (for example, chromosomes).

    Parameters
    ----------
    infile : path
        input file path or sys.stdin.buffer
    formats : list
        list of formats to pack into .cz file.
    sep :str
    usecols :list
        columns index in input file to be packed into .cz file.
    dim_cols : list
        dimensions column index, default is [0]
    chunksize :int

    Returns
    -------
    a generator, each element is a tuple, first element of tuple is a dataframe
     (header names is usecols),second element is the dims,
    in each dataframe, the dim are the same.

    """
    # cdef int i, N
    f = open1(infile)
    data, i, pre_dims = [], 0, None
    dtfuncs = get_dtfuncs(formats)
    line = f.readline()
    while line:
        values = line.rstrip('\n').split(sep)
        dims = [values[i] for i in dim_cols]
        if dims != pre_dims:  # a new dim (chrom), for example, chr1 -> chr2
            if len(data) > 0:  # write rest data of chr1
                yield pd.DataFrame(data,columns=usecols), pre_dims
                data, i = [], 0
            pre_dims = dims
        if i >= chunksize:  # dims are the same, but reach chunksize
            yield pd.DataFrame(data,columns=usecols), pre_dims
            data, i = [], 0
        data.append([func(v) for v, func in zip([values[i] for i in usecols],
                                                dtfuncs)])
        line = f.readline()
        i += 1
    f.close()
    if len(data) > 0:
        yield pd.DataFrame(data, columns=usecols), pre_dims


# ==========================================================
def _input_parser(infile, formats, sep='\t', usecols=[1, 4, 5], dim_cols=[0],
                  chunksize=5000):
    if hasattr(infile, 'readline'):
        yield from _text_input_parser(infile, formats, sep, usecols, dim_cols, chunksize)
    elif infile.endswith('.gz'):
        yield from _gz_input_parser(infile, formats, sep, usecols, dim_cols, chunksize)
    else:
        yield from _text_input_parser(infile, formats, sep, usecols, dim_cols, chunksize)


# ==========================================================
def _isCG(record):
    return record[2][:2] == b'CG'


def _isForwardCG(record):
    return record[2][:2] == b'CG' and record[1] == b'+'
# ==========================================================
def _isCH(record):
    return not _isCG(record)


# ==========================================================
class Reader:
    def __init__(self, Input, mode="rb", fileobj=None, max_cache=100):
        r"""Initialize the class for reading a BGZF file.
		You would typically use the top level ``bgzf.open(...)`` function
		which will call this class internally. Direct use is discouraged.
		Either the ``filename`` (string) or ``fileobj`` (input file object in
		binary mode) arguments must be supplied, but not both.
		Argument ``mode`` controls if the data will be returned as strings in
		text mode ("rt", "tr", or default "r"), or bytes binary mode ("rb"
		or "br"). The argument name matches the built-in ``open(...)`` and
		standard library ``gzip.open(...)`` function.
		If text mode is requested, in order to avoid multi-byte characters,
		this is hard coded to use the "latin1" encoding, and "\r" and "\n"
		are passed as is (without implementing universal new line mode). There
		is no ``encoding`` argument.

		If your data is in UTF-8 or any other incompatible encoding, you must
		use binary mode, and decode the appropriate fragments yourself.
		Argument ``max_cache`` controls the maximum number of BGZF blocks to
		cache in memory. Each can be up to 64kb thus the default of 100 blocks
		could take up to 6MB of RAM. This is important for efficient random
		access, a small value is fine for reading the file in one pass.
		"""
        if max_cache < 1:
            raise ValueError("Use max_cache with a minimum of 1")
        # Must open the BGZF file in binary mode, but we may want to
        # treat the contents as either text or binary (unicode or
        # bytes under Python 3)
        if Input and fileobj:
            raise ValueError("Supply either filename or fileobj, not both")
        # Want to reject output modes like w, a, x, +
        if mode.lower() not in ("r", "tr", "rt", "rb", "br"):
            raise ValueError(
                "Must use a read mode like 'r' (default), 'rt', or 'rb' for binary"
            )
        # If an open file was passed, make sure it was opened in binary mode.
        if fileobj:
            if fileobj.read(0) != b"":
                raise ValueError("fileobj not opened in binary mode")
            handle = fileobj
        else:
            Input=os.path.abspath(os.path.expanduser(Input))
            handle = _open(Input, "rb")
        self._handle = handle
        self.Input=Input
        self.max_cache = max_cache
        self._block_start_offset = None
        self._block_raw_length = None
        self.read_header()

    def read_header(self):
        # header: magic (5 bytes, 5s)  +
        # format_len (1byte, B) + format (format_len bytes, s)
        self.header = {}
        f = self._handle
        magic = struct.unpack("<5s", f.read(5))[0]
        if magic != _bcz_magic:
            raise ValueError("Not a right format?")
        self.header['magic'] = magic
        self.header['version'] = struct.unpack("<f", f.read(4))[0]
        total_size = struct.unpack("<Q", f.read(8))[0]
        if total_size == 0:
            raise ValueError("File not completed !")
        self.header['total_size'] = total_size
        n = struct.unpack("<H", f.read(2))[0] #message_len
        self.header['message'] = struct.unpack(f"<{n}s", f.read(n))[0].decode()
        format_len = struct.unpack("<B", f.read(1))[0]
        Formats = []
        for i in range(format_len):
            n = struct.unpack("<B", f.read(1))[0]
            format = struct.unpack(f"<{n}s", f.read(n))[0].decode()
            Formats.append(format)
        self.header['Formats'] = Formats
        Columns = []
        for i in range(format_len):
            n = struct.unpack("<B", f.read(1))[0]
            name = struct.unpack(f"<{n}s", f.read(n))[0].decode()
            Columns.append(name)
        self.header['Columns'] = Columns
        assert len(Formats) == len(Columns)
        Dimensions = []
        n_dims = struct.unpack("<B", f.read(1))[0]
        for i in range(n_dims):
            n = struct.unpack("<B", f.read(1))[0]
            dim = struct.unpack(f"<{n}s", f.read(n))[0].decode()
            Dimensions.append(dim)
        self.header['Dimensions'] = Dimensions
        self.header['header_size'] = f.tell()  # end of header, begin of 1st chunk
        self.fmts = ''.join(Formats)
        self._unit_size = struct.calcsize(self.fmts)
        self.chunk_info = self.summary_chunks(printout=False)

    def print_header(self):
        for k in self.header:
            print(k," : ",self.header[k])

    def _load_chunk(self, start_offset=None, jump=True):
        if start_offset is None:  # this is another chunk, not the 1st one.
            start_offset = self._chunk_end_offset
        if start_offset >= self.header['total_size']:
            return False
        self._handle.seek(start_offset)
        self._chunk_start_offset = start_offset  # real offset on disk.
        magic = self._handle.read(2)
        if magic != _chunk_magic:
            return False
        self._chunk_size = struct.unpack('<Q', self._handle.read(8))[0]
        # load chunk tail, jump all blocks
        self._handle.seek(self._chunk_start_offset + self._chunk_size)
        self._chunk_data_len = struct.unpack("<Q", self._handle.read(8))[0]
        self._chunk_nblocks = struct.unpack("<Q", self._handle.read(8))[0]
        self._chunk_block_1st_record_virtual_offsets = []
        if jump:  # no need to load _chunk_block_1st_record_virtual_offsets
            self._handle.seek(self._chunk_nblocks * 8, 1)
        # _chunk_tail_offset = end position of this chunk
        # = start position of next chunk.
        else:  # query,need to load block 1st record offset
            # read block_offset
            for i in range(self._chunk_nblocks):
                block_offset = struct.unpack("<Q", self._handle.read(8))[0]
                self._chunk_block_1st_record_virtual_offsets.append(block_offset)
        dims=[]
        for t in self.header['Dimensions']:
            n = struct.unpack("<B", self._handle.read(1))[0]
            dim = struct.unpack(f"<{n}s", self._handle.read(n))[0].decode()
            dims.append(dim)
        self._chunk_dims = tuple(dims)
        self._chunk_end_offset = self._handle.tell()
        return True

    def get_chunks(self):
        r = self._load_chunk(self.header['header_size'], jump=False)
        while r:
            yield [self._chunk_start_offset, self._chunk_size,
                   self._chunk_dims, self._chunk_data_len, self._chunk_end_offset,
                   self._chunk_nblocks, self._chunk_block_1st_record_virtual_offsets
                   ]
            r = self._load_chunk(jump=False)

    def summary_chunks(self, printout=True):
        r = self._load_chunk(self.header['header_size'], jump=True)
        header = ['chunk_start_offset', 'chunk_size', 'chunk_dims',
                  'chunk_tail_offset', 'chunk_nblocks', 'chunk_nrows']
        R = []
        while r:
            nrow = int(self._chunk_data_len / self._unit_size)
            record = [self._chunk_start_offset, self._chunk_size,
                      self._chunk_dims, self._chunk_end_offset,
                      self._chunk_nblocks, nrow]
            R.append(record)
            r = self._load_chunk(jump=True)
        chunk_info = pd.DataFrame(R, columns=header)
        if chunk_info.chunk_dims.value_counts().max() > 1:
            raise ValueError("Duplicated chunk dimensions detected,"
                             "Would cause conflict for querying, please check.")
        for i, dimension in enumerate(self.header['Dimensions']):
            chunk_info.insert(i, dimension, chunk_info.chunk_dims.apply(
                lambda x: x[i]
            ))
        chunk_info.set_index('chunk_dims', inplace=True)
        self.dim2chunk_start = chunk_info.chunk_start_offset.to_dict()
        if printout:
            sys.stdout.write('\t'.join(chunk_info.columns.tolist()) + '\n')
            for i, row in chunk_info.iterrows():
                try:
                    sys.stdout.write('\t'.join([str(v) for v in row.tolist()]) + '\n')
                except:
                    sys.stdout.close()
                    self.close()
                    return
            sys.stdout.close()
        else:
            return chunk_info

    def summary_blocks(self, printout=True):
        r = self._load_chunk(self.header['header_size'], jump=True)
        header = ['chunk_dims'] + ['block_start_offset', 'block_size',
                                   'block_data_start', 'block_data_len']
        if printout:
            sys.stdout.write('\t'.join(header) + '\n')
        else:
            R = []
        while r:
            self._handle.seek(self._chunk_start_offset + 10)
            chunk_info = [self._chunk_dims]
            for block in SummaryBczBlocks(self._handle):
                block = chunk_info + list(block)
                if printout:
                    try:
                        sys.stdout.write('\t'.join([str(v) for v in block]) + '\n')
                    except:
                        sys.stdout.close()
                        self.close()
                        return
                else:
                    R.append(block)
            r = self._load_chunk(jump=True)
        if printout:
            sys.stdout.close()
        else:
            df = pd.DataFrame(R, columns=header)
            return df

    def chunk2df(self, dims,reformat=False,chunksize=None):
        r = self._load_chunk(self.dim2chunk_start[dims], jump=False)
        self._cached_data = b''
        self._load_block(start_offset=self._chunk_start_offset + 10)  #
        R=[]
        i=0
        while self._block_raw_length > 0:
            # deal with such case: unit_size is 10, but data(_buffer) size is 18,
            self._cached_data += self._buffer
            end_index = len(self._cached_data) - (len(self._cached_data) % self._unit_size)
            for result in struct.iter_unpack(f"<{self.fmts}", self._cached_data[:end_index]):
                R.append(result)
            if not chunksize is None and i>= chunksize:
                df=pd.DataFrame(R,columns=self.header['Columns'])
                R=[]
                i=0
                yield df
            self._cached_data = self._cached_data[end_index:]
            self._load_block()
            i+=1
        if chunksize is None:
            df=pd.DataFrame(R,columns=self.header['Columns'])
            if not reformat:
                return df
            for col,format in zip(df.columns.tolist(),self.header['Formats']):
                # print(col,format)
                if format[-1] in ['c','s']:
                    df[col]=df[col].apply(lambda x:str(x,'utf-8'))
            return df
        else:
            if len(R) > 0:
                df = pd.DataFrame(R, columns=self.header['Columns'])
                yield df

    def _load_block(self, start_offset=None):
        if start_offset is None:
            # If the file is being read sequentially, then _handle.tell()
            # should be pointing at the start of the next block.
            # However, if seek has been used, we can't assume that.
            start_offset = self._block_start_offset + self._block_raw_length
        elif start_offset == self._block_start_offset:
            self._within_block_offset = 0
            return
        # Now load the block
        self._handle.seek(start_offset)
        self._block_start_offset = start_offset
        try:
            block_size, self._buffer = _load_bcz_block(self._handle, True)
        except StopIteration:  # EOF
            block_size = 0
            self._buffer = b""
        self._within_block_offset = 0
        self._block_raw_length = block_size

    def _byte2str(self, values):
        return [str(v, 'utf-8') if f[-1] in ['s', 'c'] else str(v)
                for v, f in zip(values, self.header['Formats'])]

    def _byte2real(self, values):
        return [str(v, 'utf-8') if f[-1] in ['s', 'c'] else v
                for v, f in zip(values, self.header['Formats'])]

    def _empty_generator(self):
        while True:
            yield []

    def view(self, show_dim=None, header=True, Dimension=None,
             reference=None):
        """
		View .cz file.

		Parameters
		----------
		show_dim : str
			index of dims given to writer.write_chunk, separated by comma,
			default is None, dims[show_dim] will be shown in each row (such as sampleID
			and chrom)
		header : bool
			whether to print the header.
        Dimension: None, bool, list or file path
            If None (default): use the default chunk order in .cz file; \
            If list: use this Dimension (dims) as order and print records.\
            If path: use the first len(Dimension) columns as dim order, there should be\
                no header in file path and use \t as separator.\
            If Dimension is dictionary,such as -D "{'sampleID':'cell1','chrom':'chr1'}",\
            will filter chunk using sampleID=cell1 and chrom=chr1.

		Returns
		-------

		"""
        if isinstance(show_dim, int):
            show_dim = [show_dim]
        elif isinstance(show_dim, str):
            show_dim = [int(i) for i in show_dim.split(',')]

        if Dimension is None or isinstance(Dimension, dict):
            chunk_info = self.chunk_info.copy()
            if isinstance(Dimension, dict):
                selected_dim = Dimension.copy()
                for d, v in selected_dim.items():
                    chunk_info = chunk_info.loc[chunk_info[d] == v]
            Dimension = chunk_info.index.tolist()
        else:
            # equal to query chromosome if self.header['dimensions'][0]==chrom
            if isinstance(Dimension, str):
                order_path = os.path.abspath(os.path.expanduser(Dimension))
                if os.path.exists(order_path):
                    dim = pd.read_csv(order_path, sep='\t', header=None,
                                      usecols=show_dim)[show_dim].apply(lambda x: tuple(x.tolist()),
                                                                        axis=1).tolist()
                else:
                    Dimension = Dimension.split(',')
            if isinstance(Dimension, (list, tuple)) and isinstance(Dimension[0], str):
                Dimension = [tuple([o]) for o in Dimension]
            if not isinstance(Dimension, (list, tuple, np.ndarray)):  # dim is a list
                raise ValueError("input of dim_order is not corrected !")

        if not reference is None:
            reference = os.path.abspath(os.path.expanduser(reference))
            ref_reader = Reader(reference)
        if not show_dim is None:
            header_columns = [self.header['Dimensions'][t] for t in show_dim]
            dim_header = "\t".join(header_columns) + '\t'
        else:
            dim_header = ''
        if header:  # show header
            if not reference is None:
                columns = ref_reader.header['Columns'] + self.header['Columns']
            else:
                columns = self.header['Columns']
            line = "\t".join(columns)
            sys.stdout.write(dim_header + line + '\n')

        for d in Dimension:  # Dimension is a list of tuple ([(d1,d2),(d1,d2)])
            if d not in self.dim2chunk_start:
                continue
            if not show_dim is None:
                dim_stdout = "\t".join([d[t] for t in show_dim]) + '\t'
            else:
                dim_stdout = ''
            records = self.fetch(d)
            if not reference is None:
                ref_records = ref_reader.fetch(d)
            else:
                ref_records = self._empty_generator()
            try:
                for record, ref_record in zip(records, ref_records):
                    line = '\t'.join([str(i) for i in ref_record + record])
                    try:
                        sys.stdout.write(dim_stdout + line + '\n')
                    except:
                        sys.stdout.close()
                        if not reference is None:
                            ref_reader.close()
                        self.close()
                        return
            except:
                raise ValueError(f"reference {reference} not matched.")
        sys.stdout.close()
        self.close()
        if not reference is None:
            ref_reader.close()

    def category_ssi(self, output=None, formats=['I'], columns=['ID'],
                     dimensions=['chrom'], match_func=_isForwardCG,
                     chunksize=2000):
        if output is None:
            output = self.Input + '.' + match_func.__name__ + '.ssi'
        else:
            output = os.path.abspath(os.path.expanduser(output))
        writer = Writer(output, Formats=formats, Columns=columns,
                        Dimensions=dimensions, fileobj=None,
                        message=os.path.basename(self.Input))
        data = b''
        for dim in self.dim2chunk_start:
            print(dim, '\t' * 4, end='\r')
            for i, record in enumerate(self.__fetch__(dim)):
                if match_func(record):
                    data += struct.pack(f"<{writer.fmts}", i + 1)
                if (i % chunksize) == 0 and len(data) > 0:
                    writer.write_chunk(data, dim)
                    data = b''
            if len(data) > 0:
                writer.write_chunk(data, dim)
                data = b''
        writer.close()

    def get_ids_from_ssi(self, dim):
        if len(self.header['Columns']) == 1:
            s, e = 0, 1  # only one columns, ID
            return np.array([record[0] for record in self.__fetch__(dim, s=s, e=e)])
        else:  # two columns: ID_start, ID_end
            s, e = 0, 2
            return np.array([record for record in self.__fetch__(dim, s=s, e=e)])

    def _getRecordsByIds(self, dim=None, IDs=None):
        """

        Parameters
        ----------
        dim : tuple
        reference : path
        IDs : np.array

        Returns
        -------

        """
        self._load_chunk(self.dim2chunk_start[dim], jump=False)
        block_index = ((IDs - 1) * self._unit_size) // (_BLOCK_MAX_LEN)
        block_start_offsets = np.array([self._chunk_block_1st_record_virtual_offsets[
                                            idx] >> 16 for idx in block_index])
        within_block_offsets = ((IDs - 1) * self._unit_size) % (_BLOCK_MAX_LEN)
        block_start_offset_tmp = None
        for block_start_offset, within_block_offset in zip(
            block_start_offsets, within_block_offsets):
            if block_start_offset != block_start_offset_tmp:
                self._load_block(block_start_offset)
                block_start_offset_tmp = block_start_offset
            self._within_block_offset = within_block_offset
            yield self.read(self._unit_size)

    def getRecordsByIds(self, dim=None, reference=None, IDs=None):
        if reference is None:
            for record in self._getRecordsByIds(dim, IDs):
                yield self._byte2real(struct.unpack(
                    f"<{self.fmts}", record
                ))
        else:
            ref_reader = Reader(os.path.abspath(os.path.expanduser(reference)))
            ref_records = ref_reader._getRecordsByIds(dim=dim, IDs=IDs)
            records = self._getRecordsByIds(dim=dim, IDs=IDs)
            for ref_record, record in zip(ref_records, records):
                yield ref_reader._byte2real(struct.unpack(
                    f"<{ref_reader.fmts}", ref_record
                )) + self._byte2real(struct.unpack(
                    f"<{self.fmts}", record
                ))
            ref_reader.close()

    def _getRecordsByIdRegions(self, dim=None, IDs=None):
        """
        Get .cz content for a given dim and IDs
        Parameters
        ----------
        dim : tuple
        reference : path
        IDs : list or np.ndarray
            every element of IDs is a list or tuple with length=2, id_start and id_end

        Returns
        -------
        generator
        """
        self._load_chunk(self.dim2chunk_start[dim], jump=False)
        block_start_offsets_tmp = None
        for id_start, id_end in IDs:
            block_index = ((id_start - 1) * self._unit_size) // (_BLOCK_MAX_LEN)
            block_start_offset = self._chunk_block_1st_record_virtual_offsets[
                                     block_index] >> 16
            if block_start_offset != block_start_offsets_tmp:
                self._load_block(block_start_offset)
                block_start_offsets_tmp = block_start_offset
            self._within_block_offset = ((id_start - 1) * self._unit_size
                                         ) % _BLOCK_MAX_LEN
            yield [self.read(self._unit_size) for i in range(id_start, id_end + 1)]

    def getRecordsByIdRegions(self, dim=None, reference=None, IDs=None):
        """
        Get .cz content for a given dim and IDs
        Parameters
        ----------
        dim : tuple
        reference : path
        IDs : list or np.ndarray
            every element of IDs is a list or tuple with length=2, id_start and id_end

        Returns
        -------
        generator
        """
        self._load_chunk(self.dim2chunk_start[dim], jump=False)
        if reference is None:
            for records in self._getRecordsByIdRegions(dim, IDs):
                yield np.array([self._byte2real(struct.unpack(
                    f"<{self.fmts}", record)) for record in records])
        else:
            ref_reader = Reader(os.path.abspath(os.path.expanduser(reference)))
            ref_records = ref_reader._getRecordsByIdRegions(dim=dim, IDs=IDs)
            records = self._getRecordsByIdRegions(dim=dim, IDs=IDs)
            for ref_records, records in zip(ref_records, records):
                ref_record = np.array([ref_reader._byte2real(struct.unpack(
                    f"<{ref_reader.fmts}", record)) for record in ref_records])
                record = np.array([self._byte2real(struct.unpack(
                    f"<{self.fmts}", record)) for record in records])
                yield np.hstack((ref_record, record))
            ref_reader.close()

    def subset(self, dim, ssi=None, IDs=None, reference=None, printout=True):
        if isinstance(dim, str):
            dim = tuple([dim])
        if ssi is None and IDs is None:
            raise ValueError("Please provide either ssi or IDs")
        if not ssi is None:
            ssi_reader = Reader(ssi)
            IDs = ssi_reader.get_ids_from_ssi(dim)
            ssi_reader.close()
        if not reference is None:
            ref_reader = Reader(os.path.abspath(os.path.expanduser(reference)))
            ref_header = ref_reader.header['Columns']
            ref_reader.close()
        else:
            ref_header = []
        header = self.header['Dimensions'] + ref_header + self.header['Columns']
        if len(IDs.shape) == 1:
            if printout:
                sys.stdout.write('\t'.join(header) + '\n')
                try:
                    for record in self.getRecordsByIds(dim, reference, IDs):
                        sys.stdout.write('\t'.join(list(dim) + list(map(str, record))) + '\n')
                        # print(list(dim)+list(map(str,record)))
                except:
                    sys.stdout.close()
                    self.close()
                    return
                sys.stdout.close()
            else:  # each element is one record
                yield from self.getRecordsByIds(dim, reference, IDs)
        else:
            if printout:
                sys.stdout.write('\t'.join(header) + '\n')
                try:
                    for data in self.getRecordsByIdRegions(dim, reference, IDs):
                        for row in data:
                            sys.stdout.write('\t'.join(list(dim) + list(map(str, row))) + '\n')
                except:
                    sys.stdout.close()
                    self.close()
                    return
                sys.stdout.close()
            else:  # each element is a data array (multiple records)
                yield from self.getRecordsByIdRegions(dim, reference, IDs)

    def __fetch_deprecated__(self, dims, s=None, e=None):
        """
        Generator for a given dims

        Parameters
        ----------
        dims : tuple
            element length of dims should match the Dimensions in header['Dimensions']

        Returns
        -------

        """
        s = 0 if s is None else s  # 0-based
        e = len(self.header['Columns']) if e is None else e  # 1-based
        r = self._load_chunk(self.dim2chunk_start[dims], jump=True)
        self._cached_data = b''
        # unit_nblock = int(self._unit_size / (math.gcd(self._unit_size, _BLOCK_MAX_LEN)))
        self._load_block(start_offset=self._chunk_start_offset + 10)  #
        while self._block_raw_length > 0:
            self._cached_data += self._buffer
            end_index = len(self._cached_data) - (len(self._cached_data) % self._unit_size)
            for result in struct.iter_unpack(f"<{self.fmts}", self._cached_data[:end_index]):
                yield result[s:e]  # a tuple
                # print(result)
            self._cached_data = self._cached_data[end_index:]
            self._load_block()

    def __fetch__(self, dims, s=None, e=None):
        """
        Generator for a given dims

        Parameters
        ----------
        dims : tuple
            element length of dims should match the Dimensions in header['Dimensions']

        Returns
        -------

        """
        s = 0 if s is None else s  # 0-based
        e = len(self.header['Columns']) if e is None else e  # 1-based
        r = self._load_chunk(self.dim2chunk_start[dims], jump=True)
        unit_nblock = int(self._unit_size / (math.gcd(self._unit_size, _BLOCK_MAX_LEN)))
        self._load_block(start_offset=self._chunk_start_offset + 10)  #
        i, _cached_data = 1, self._buffer
        while self._block_raw_length > 0:
            self._load_block()
            if i < unit_nblock:
                _cached_data += self._buffer
                i += 1
            else:
                for result in struct.iter_unpack(f"<{self.fmts}", _cached_data):
                    yield result[s:e]  # a tuple
                i, _cached_data = 1, self._buffer
        if len(_cached_data) > 0:
            for result in struct.iter_unpack(f"<{self.fmts}", _cached_data):
                yield result[s:e]  # a tuple

    def fetch(self, dims):
        for record in self.__fetch__(dims):
            yield self._byte2real(record)  # all columns of each row

    def batch_fetch(self, dims, chunksize=5000):
        i = 0
        data = []
        for row in self.fetch(dims):
            data.append(row)
            if i >= chunksize:
                yield data
                data, i = [], 0
            i += 1
        if len(data) > 0:
            yield data

    def fetchByStartID(self, dims, n=None):  # n is 1-based, n >=1
        if isinstance(dims, str):
            dims = tuple([dims])
        r = self._load_chunk(self.dim2chunk_start[dims], jump=False)
        if not n is None:  # seek to the position of the n rows
            block_index = ((n - 1) * self._unit_size) // (_BLOCK_MAX_LEN)
            first_record_vos = self._chunk_block_1st_record_virtual_offsets[block_index]
            within_block_offset = ((n - 1) * self._unit_size) % (_BLOCK_MAX_LEN)
            block_start_offset = (first_record_vos >> 16)
            virtual_start_offset = (block_start_offset << 16) | within_block_offset
            self.seek(virtual_start_offset)  # load_block is inside sek
        while True:
            # yield self._read_1record()
            yield struct.unpack(f"<{self.fmts}", self.read(self._unit_size))

    def _read_1record(self):
        tmp = self._buffer[
              self._within_block_offset:self._within_block_offset + self._unit_size]
        self._within_block_offset += self._unit_size
        return struct.unpack(f"<{self.fmts}", tmp)

    def _query_regions(self, regions, s, e):
        """
        query given regions, start and end index to the columns of data columns.

        Parameters
        ----------
        regions : list
            list of list, first element is a tuple (dims), second is start, third is end,.
        s : int
        e : int

        Returns
        -------
        Return a generator, the first element of generator is start offset of the
        1st record.
        """
        dim_tmp = None
        for dim, start, end in regions:
            if dim != dim_tmp:
                # start_block_index_tmp = 0
                start_block_index = 0
                dim_tmp = dim
                r = self._load_chunk(self.dim2chunk_start[dim], jump=False)
                # find the closest block to a given start position
                block_1st_starts = np.array([
                    self._seek_and_read_1record(offset)[s]
                    for offset in self._chunk_block_1st_record_virtual_offsets
                ])
            # for idx in range(start_block_index_tmp, self._chunk_nblocks-1):
            #     if block_1st_starts[idx + 1] > start:
            #         start_block_index = idx
            #         break
            #     start_block_index = self._chunk_nblocks - 1
            while start_block_index < self._chunk_nblocks - 1:
                if block_1st_starts[start_block_index + 1] > start:
                    break
                start_block_index += 1
            # start_block_index=np.argmax(block_1st_starts > start)-1
            virtual_offset = self._chunk_block_1st_record_virtual_offsets[start_block_index]
            self.seek(virtual_offset)  # seek to the target block, self._buffer
            block_start_offset = self._block_start_offset
            record = struct.unpack(f"<{self.fmts}", self.read(self._unit_size))
            while record[s] < start:
                # record = self._read_1record()
                record = struct.unpack(f"<{self.fmts}", self.read(self._unit_size))
            # here we got the 1st record, return the id (row number) of 1st record.
            # size of each block is 65535 (2**16-1=_BLOCK_MAX_LEN)
            # primary_id = int((_BLOCK_MAX_LEN * block_index +
            #                   self._within_block_offset) / self._unit_size)
            if self._block_start_offset > block_start_offset:
                start_block_index += 1
            primary_id = int((_BLOCK_MAX_LEN * start_block_index +
                              self._within_block_offset) / self._unit_size)
            # primary_id is the current record (record[s] == start),
            # with_block_offset is the end of current record
            # 1-based, primary_id >=1, cause _within_block_offset >= self._unit_size
            yield "primary_id_&_dim:", primary_id, dim
            while record[e] <= end:
                yield dim, record
                record = struct.unpack(f"<{self.fmts}", self.read(self._unit_size))
            start_block_index_tmp = start_block_index

    def pos2id(self, dim, positions, col_to_query=0):  # return IDs (primary_id)
        self._load_chunk(self.dim2chunk_start[dim], jump=False)
        block_1st_starts = [
            self._seek_and_read_1record(offset)[col_to_query]
            for offset in self._chunk_block_1st_record_virtual_offsets
        ]
        start_block_index_tmp = 0
        R = []
        for pos in positions:  # positions is sorted, [[0,1],[3,4],[7,8]]
            for idx in range(start_block_index_tmp, self._chunk_nblocks - 1):
                if block_1st_starts[idx + 1] > pos:
                    start_block_index = idx
                    break
                start_block_index = self._chunk_nblocks - 1
            virtual_offset = self._chunk_block_1st_record_virtual_offsets[start_block_index]
            self.seek(virtual_offset)  # seek to the target block, self._buffer
            block_start_offset = self._block_start_offset
            record = struct.unpack(f"<{self.fmts}", self.read(self._unit_size))
            while record[col_to_query] < pos:
                record = struct.unpack(f"<{self.fmts}", self.read(self._unit_size))
            if self._block_start_offset > block_start_offset:
                start_block_index += 1
            primary_id = int((_BLOCK_MAX_LEN * start_block_index +
                              self._within_block_offset) / self._unit_size)
            R.append(primary_id)
            start_block_index_tmp = start_block_index
        return R

    def query(self, Dimension=None, start=None, end=None, Regions=None,
              query_col=[0], reference=None, printout=True):
        """
        query .cz file by Dimension, start and end, if regions provided, Dimension, start and
        end should be None, regions should be a list, each element of regions
        is a list, for example, regions=[[('cell1','chr1'),1,10],
        [('cell10','chr22'),100,200]],and so on.

        Parameters
        ----------
        Dimension : str
            dimension, such as chr1 (if Dimensions 'chrom' is inlcuded in
            header['Dimensions']), or sample1 (if something like 'sampleID' is
            included in header['Dimensions'] and chunk contains such dimension)
        start : int
            start position, if None, the entire Dimension would be returned.
        end : int
            end position
        regions : list or file path
            regions=[
                [Dimension,start,end], [Dimension,start,end]
            ],
            Dimension is a tuple, for example, Dimension=('chr1');

            if regions is a file path (separated by tab and no header),
            the first len(header['Dimensions']) columns would be
            used as Dimension, and next two columns would be used as start and end.
        query_col: list
            index of columns (header['Columns']) to be queried,for example,
            if header['Columns']=['pos','mv','cov'], then pos is what we want to
            query, so query_col should be [0], but if
            header['Columns']=['start','end','peak'], then start and end are what
            we would like to query, query_col should be [0,1]

        Returns
        -------

        """
        if Dimension is None and Regions is None:
            raise ValueError("Please provide either Dimension,start,end or regions")
        if (not Dimension is None) and (not Regions is None):
            raise ValueError("Please query by Dimension,start,end or by regions, "
                             "not both !")
        if not Dimension is None:
            assert not start is None
            assert not end is None
            if isinstance(Dimension, str):
                Dimension = tuple([Dimension])
                Regions = [[Dimension, start, end]]
            elif isinstance(Dimension, dict):
                chunk_info = self.chunk_info.copy()
                selected_dim = Dimension.copy()
                for d, v in selected_dim.items():
                    chunk_info = chunk_info.loc[chunk_info[d] == v]
                Dimension = chunk_info.index.tolist()
                Regions = [[dim, start, end] for dim in Dimension]
        else:
            if isinstance(Regions, str):
                region_path = os.path.abspath(os.path.expanduser(Regions))
                if os.path.exists(region_path):
                    n_dim = self.header['Dimensions']
                    usecols = list(range(n_dim + 2))
                    df = pd.read_csv(region_path, sep='\t', header=None,
                                     usecols=usecols)
                    Regions = df.apply(lambda x: [
                        tuple(x[:n_dim]),
                        x[n_dim],
                        x[n_dim + 1]
                    ], axis=1).tolist()
                else:
                    raise ValueError(f"path {region_path} not existed.")
            else:  # check format of regions
                if not isinstance(Regions[0][0], tuple):
                    raise ValueError("The first element of first region is not a tuple")

        if len(query_col) == 1:
            s = e = query_col[0]
        elif len(query_col) == 2:
            s, e = query_col
        else:
            raise ValueError("length of query_col can not be greater then 2.")

        if reference is None:  # query position in this .cz file.
            for i in set([s, e]):
                if self.header['Formats'][i][-1] in ['s', 'c']:
                    raise ValueError("The query_col of Columns is not int or float",
                                     "Can not perform querying without reference!",
                                     f"Columns: {self.header['Columns']}",
                                     f"Formats: {self.header['Formats']}")
            header = self.header['Dimensions'] + self.header['Columns']
            if printout:
                sys.stdout.write('\t'.join(header) + '\n')
                records = self._query_regions(Regions, s, e)
                record = records.__next__()
                while True:
                    try:
                        record = records.__next__()
                    except:
                        break
                    try:
                        dims, row = record
                        try:
                            sys.stdout.write('\t'.join([str(d) for d in dims] +
                                                       self._byte2str(row)) + '\n')
                        except:
                            sys.stdout.close()
                            self.close()
                            return
                    except ValueError:  # len(record) ==3: #"primary_id_&_dim:", primary_id, dim
                        flag, primary_id, dim = record
                sys.stdout.close()
                self.close()
            else:
                records = self._query_regions(Regions, s, e)
                record = records.__next__()
                while True:
                    try:
                        record = records.__next__()
                    except:
                        break
                    try:
                        dims, row = record
                        yield list(dims) + list(self._byte2real(row))
                    except ValueError:
                        flag, primary_id, dim = record
        else:
            reference = os.path.abspath(os.path.expanduser(reference))
            ref_reader = Reader(reference)
            header = (self.header['Dimensions'] + ref_reader.header['Columns']
                      + self.header['Columns'])
            if printout:
                sys.stdout.write('\t'.join(header) + '\n')
                # query reference first
                ref_records = ref_reader._query_regions(Regions, s, e)
                flag, primary_id, dim = ref_records.__next__()  # 1st should contain primary_id
                records = self.fetchByStartID(dim, n=primary_id)
                while True:
                    try:
                        ref_record = next(ref_records)
                        record = self._byte2str(next(records))
                    except:
                        break
                    try:
                        dims, row = ref_record
                        rows = ref_reader._byte2str(row) + record
                        try:
                            sys.stdout.write('\t'.join([str(d) for d in dims] + rows) + '\n')
                        except:
                            sys.stdout.close()
                            self.close()
                            return
                    except ValueError:  # len(record) ==3: #"primary_id_&_dim:", primary_id, dim
                        flag, primary_id, dim = ref_record  # next block or next dim
                        records = self.fetchByStartID(dim, n=primary_id)
                sys.stdout.close()
                self.close()
            else:
                ref_records = ref_reader._query_regions(Regions, s, e)
                ref_record = ref_records.__next__()  # 1st should contain primary_id
                flag, primary_id, dim = ref_record
                records = self.fetchByStartID(dim, n=primary_id)
                while True:
                    try:
                        ref_record = ref_records.__next__()
                        record = self._byte2str(records.__next__())
                    except:
                        break
                    try:
                        dims, row = ref_record
                        rows = ref_reader._byte2str(row) + record
                        yield list(dims) + rows
                    except ValueError:  # len(record) ==3: #"primary_id_&_dim:", primary_id, dim
                        flag, primary_id, dim = ref_record  # next block or next dim
                        records = self.fetchByStartID(dim, n=primary_id)
            ref_reader.close()

    def _seek_and_read_1record(self,virtual_offset):
        self.seek(virtual_offset)
        return struct.unpack(f"<{self.fmts}", self.read(self._unit_size))

    def tell(self):
        """Return a 64-bit unsigned BGZF virtual offset."""
        if 0 < self._within_block_offset and self._within_block_offset == len(
            self._buffer):
            # Special case where we're right at the end of a (non empty) block.
            # For non-maximal blocks could give two possible virtual offsets,
            # but for a maximal block can't use _BLOCK_MAX_LEN as the within block
            # offset. Therefore for consistency, use the next block and a
            # within block offset of zero.
            return (self._block_start_offset + self._block_raw_length) << 16
        else:
            return (self._block_start_offset << 16) | self._within_block_offset

    def seek(self, virtual_offset):
        """Seek to a 64-bit unsigned BGZF virtual offset."""
        # Do this inline to avoid a function call,
        # start_offset, within_block = split_virtual_offset(virtual_offset)
        start_offset = virtual_offset >> 16
        within_block = virtual_offset ^ (start_offset << 16)
        if start_offset != self._block_start_offset:
            # Don't need to load the block if already there
            self._load_block(start_offset)
            if start_offset != self._block_start_offset:
                raise ValueError("start_offset not loaded correctly")
        if within_block > len(self._buffer):
            if not (within_block == 0 and len(self._buffer) == 0):
                raise ValueError(
                    "Within offset %i but block size only %i"
                    % (within_block, len(self._buffer))
                )
        self._within_block_offset = within_block
        return virtual_offset

    def read(self, size=1):
        """Read method for the BGZF module."""
        result = b""
        while size and self._block_raw_length:
            if self._within_block_offset + size <= len(self._buffer):
                # This may leave us right at the end of a block
                # (lazy loading, don't load the next block unless we have too)
                data = self._buffer[
                       self._within_block_offset: self._within_block_offset + size
                       ]
                self._within_block_offset += size
                # if not data:
                #     raise ValueError("Must be at least 1 byte")
                result += data
                break
            else:  # need to load the enxt block
                data = self._buffer[self._within_block_offset:]
                size -= len(data)
                self._load_block()  # will reset offsets and _within_block_offset
                result += data

        return result

    def readline(self):
        """Read a single line for the BGZF file."""
        result = b""
        while self._block_raw_length:
            i = self._buffer.find(self._newline, self._within_block_offset)
            # Three cases to consider,
            if i == -1:
                # No newline, need to read in more data
                data = self._buffer[self._within_block_offset:]
                self._load_block()  # will reset offsets
                result += data
            elif i + 1 == len(self._buffer):
                # Found new line, but right at end of block (SPECIAL)
                data = self._buffer[self._within_block_offset:]
                # Must now load the next block to ensure tell() works
                self._load_block()  # will reset offsets
                if not data:
                    raise ValueError("Must be at least 1 byte")
                result += data
                break
            else:
                # Found new line, not at end of block (easy case, no IO)
                data = self._buffer[self._within_block_offset: i + 1]
                self._within_block_offset = i + 1
                # assert data.endswith(self._newline)
                result += data
                break

        return result

    def close(self):
        """Close BGZF file."""
        self._handle.close()
        self._buffer = None
        self._block_start_offset = None
        self._buffers = None

    def seekable(self):
        """Return True indicating the BGZF supports random access."""
        return True

    def isatty(self):
        """Return True if connected to a TTY device."""
        return False

    def fileno(self):
        """Return integer file descriptor."""
        return self._handle.fileno()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()


# ==========================================================
def extract(input=None, outfile=None, ssi=None, chunksize=5000):
    ssi_reader = Reader(os.path.abspath(os.path.expanduser(ssi)))
    reader = Reader(os.path.abspath(os.path.expanduser(input)))
    writer = Writer(outfile, Formats=reader.header['Formats'],
                    Columns=reader.header['Columns'],
                    Dimensions=reader.header['Dimensions'],
                    message=ssi)
    # dtfuncs = get_dtfuncs(writer.Formats)
    for dim in reader.dim2chunk_start.keys():
        print(dim)
        IDs = ssi_reader.get_ids_from_ssi(dim)
        if len(IDs.shape) != 1:
            raise ValueError("Only support 1D ssi now!")
        records = reader._getRecordsByIds(dim, IDs)
        data, i = b'', 0
        for record in records:  # unpacked bytes
            data += record
            i += 1
            if i > chunksize:
                writer.write_chunk(data, dim)
                data, i = b'', 0
        if len(data) > 0:
            writer.write_chunk(data, dim)
    writer.close()
    reader.close()
    ssi_reader.close()


# ==========================================================
class Writer:
    def __init__(self, Output=None, mode="wb", Formats=['H', 'H'],
                 Columns=['mc', 'cov'], Dimensions=['chrom'], fileobj=None,
                 message='', level=6, verbose=0):
        """
        czip .cz writer.

        Parameters
        ----------
        Output : path or None
            if None (default) bytes stream would be written to stdout.
        mode : str
            'wb', 'ab'
        Formats : list
            format for each column, see https://docs.python.org/3/library/struct.html#format-characters for detail format.
        Columns : list
            columns names, length should be the same as Formats.
        Dimensions : list
            Dimensions to be included for each chunk, for example, if you would like
            to include sampleID and chromosomes as dimension title, then set
            Dimsnsion=['sampleID','chrom'], then give each chunk a dims,
            for example, dims for chunk1 is ['cell1','chr1'], dims for chunk2 is
            ['cell1','chr2'], ['cell2','chr1'],....['celln','chr22']...
        fileobj : object
            an openned file object
        message : str
            a message to be included in the header, a genome assemble version
            is highly recommended to be set as message.
        level : int
            compress level (default is 6)
        verbose : int
            whether to print debug information
        """
        if Output and fileobj:
            raise ValueError("Supply either Output or fileobj, not both")
        if fileobj:  # write to an existed openned file object
            if fileobj.read(0) != b"":
                raise ValueError("fileobj not opened in binary mode")
            handle = fileobj
        else:
            if not Output is None:  # write output to a file
                if "w" not in mode.lower() and "a" not in mode.lower():
                    raise ValueError(f"Must use write or append mode")
                handle = _open(Output, mode)
                self.Output=Output
            else:  # write to stdout buffer
                handle = sys.stdout.buffer
        self._handle = handle
        self._buffer = b""
        self._chunk_start_offset = None
        self._chunk_dims = None
        if isinstance(Formats, str) and ',' not in Formats:
            self.Formats = [Formats]
        elif isinstance(Formats, str):
            self.Formats = Formats.split(',')
        else:
            self.Formats = list(Formats)
        self.level = level
        if isinstance(Columns, str):
            self.Columns = Columns.split(',')
        else:
            self.Columns = list(Columns)
        if isinstance(Dimensions, str):
            self.Dimensions = Dimensions.split(',')
        else:
            self.Dimensions = list(Dimensions)
        self._magic_size = len(_bcz_magic)
        self.verbose = verbose
        self.message = message
        # if self.verbose > 0:
        #     print(self.Formats, self.Columns, self.Dimensions)
        #     print(type(self.Formats), type(self.Columns), type(self.Dimensions))
        self.write_header()

    def write_header(self):
        f = self._handle
        f.write(struct.pack(f"<{self._magic_size}s", _bcz_magic))  # 5 bytes
        f.write(struct.pack("<f", _version))  # 4 bytes, float
        f.write(struct.pack("<Q", 0))  # 8 bytes, total size, including magic.
        f.write(struct.pack("<H", len(self.message)))  # length of each format, 1 byte
        f.write(struct.pack(f"<{len(self.message)}s",
                            bytes(self.message, 'utf-8')))
        f.write(struct.pack("<B", len(self.Formats)))  # 1byte, ncols
        assert len(self.Columns) == len(self.Formats)
        for format in self.Formats:
            format_len = len(format)
            f.write(struct.pack("<B", format_len))  # length of each format, 1 byte
            f.write(struct.pack(f"<{format_len}s", bytes(format, 'utf-8')))
        for name in self.Columns:
            name_len = len(name)
            f.write(struct.pack("<B", name_len))  # length of each name, 1 byte
            f.write(struct.pack(f"<{name_len}s", bytes(name, 'utf-8')))
        self._n_dim_offset=f.tell()
        #when a new dim is added, go back here and rewrite the n_dim (1byte)
        f.write(struct.pack("<B", len(self.Dimensions)))  # 1byte
        for dim in self.Dimensions:
            dname_len = len(dim)
            f.write(struct.pack("<B", dname_len))  # length of each name, 1 byte
            f.write(struct.pack(f"<{dname_len}s", bytes(dim, 'utf-8')))
        # when multiple .cz are cat into one .cz and a new dimension is added,
        # such as sampleID, seek to this position (_header_size) and write another
        # two element: new_dname_len (B) and new_dim, then go to
        # _n_dim_offset,rewrite the n_dim (n_dim = n_dim + 1) and
        # self.Dimensions.append(new_dim)
        self._header_size=f.tell()
        self.fmts = ''.join(list(self.Formats))
        self._unit_size = struct.calcsize(self.fmts)

    def _write_block(self, block):
        # if len(block) > _BLOCK_MAX_LEN:  # 65536 = 1 << 16 = 2**16
        #     raise ValueError(f"{len(block)} Block length > {_BLOCK_MAX_LEN}")
        c = zlib.compressobj(
            self.level, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
        )
        """
		deflate_compress = zlib.compressobj(9, zlib.DEFLATED, -zlib.MAX_WBITS)
		zlib_compress = zlib.compressobj(9, zlib.DEFLATED, zlib.MAX_WBITS)
		gzip_compress = zlib.compressobj(9, zlib.DEFLATED, zlib.MAX_WBITS | 16)
		"""
        compressed = c.compress(block) + c.flush()
        del c
        # if len(compressed) > _BLOCK_MAX_LEN:
        #     raise RuntimeError(
        #         "Didn't compress enough, try less data in this block"
        #     )
        bsize = struct.pack("<H", len(compressed) + 6)
        # block size: magic (2 btyes) + block_size (2bytes) + compressed data +
        # block_data_len (2 bytes)
        data_len = len(block)
        uncompressed_length = struct.pack("<H", data_len)  # 2 bytes
        data = _block_magic + bsize + compressed + uncompressed_length
        block_start_offset = self._handle.tell()
        within_block_offset = int(
            np.ceil(self._chunk_data_len / self._unit_size) * self._unit_size
            - self._chunk_data_len)
        virtual_offset = (block_start_offset << 16) | within_block_offset
        self._block_1st_record_virtual_offsets.append(virtual_offset)
        # _block_offsets are real position on disk, not a virtual offset
        self._handle.write(data)
        # how many bytes (uncompressed) have been writen.
        self._chunk_data_len += data_len

    def _write_chunk_tail(self):
        # tail: chunk_data_len (Q,8bytes) + n_blocks (Q, 8bytes)
        # + list of block start offsets (n_blocks * 8 bytes)
        # + List of dimensions ()
        self._handle.write(struct.pack("<Q", self._chunk_data_len))
        self._handle.write(struct.pack("<Q", len(self._block_1st_record_virtual_offsets)))
        for block_offset in self._block_1st_record_virtual_offsets:
            self._handle.write(struct.pack("<Q", block_offset))
        # write list of dims
        for dim in self._chunk_dims:
            dim_len = len(dim)
            self._handle.write(struct.pack("<B", dim_len))  # length of each dim, 1 byte
            self._handle.write(struct.pack(f"<{dim_len}s", bytes(dim, 'utf-8')))

    def _chunk_finished(self):
        # current position is the end of chunk.
        cur_offset = self._handle.tell()  # before chunk tail, not a virtual offset
        # go back and write the total chunk size (tail not included)
        self._handle.seek(self._chunk_start_offset + 2)  # 2 bytes for magic
        chunk_size = cur_offset - self._chunk_start_offset
        # chunk_size including the chunk_size itself, but not including chunk tail.
        self._handle.write(struct.pack("<Q", chunk_size))
        # go back the current position
        self._handle.seek(cur_offset)  # not a virtual offset
        self._write_chunk_tail()

    def write_chunk(self, data, dims):  # dims is a list.
        """
        Write data into one chunk with dimension name (such as ['cell1','chr1'])

        Parameters
        ----------
        data : bytes
            In general data should be bytes.
        dims : list or tuple
            dimensions to be written to chunk, length of dims should be the same
            as Dimensions given to Writer.

        Returns
        -------

        """
        # assert len(dims)==len(self.Dimensions)
        # if isinstance(data, str):
        #     data = data.encode("latin-1")
        if self._chunk_dims != dims:
            # the first chunk or another new chunk
            if self._chunk_dims is not None:
                # this is another new chunk, current position is the end of chunk.
                self.flush()  # finish_chunk in flush
            # else: #the first chunk, no chunk was writen previously.
            self._chunk_dims = dims
            self._chunk_start_offset = self._handle.tell()
            self._handle.write(_chunk_magic)
            # chunk total size place holder: 0
            self._handle.write(struct.pack("<Q", 0))  # 8 bytes; including this chunk_size
            self._chunk_data_len = 0
            self._block_1st_record_virtual_offsets = []
            # if self.verbose > 0:
            #     print("Writing chunk with dims: ", self._chunk_dims)

        if len(self._buffer) + len(data) < _BLOCK_MAX_LEN:
            self._buffer += data
        else:
            self._buffer += data
            while len(self._buffer) >= _BLOCK_MAX_LEN:
                self._write_block(self._buffer[:_BLOCK_MAX_LEN])
                self._buffer = self._buffer[_BLOCK_MAX_LEN:]

    def _parse_input_no_ref(self, input_handle):
        if isinstance(input_handle, pd.DataFrame):
            # usecols and dim_cols should be in the columns of this dataframe.
            if self.chunksize is None:
                for dim, df1 in input_handle.groupby(self.dim_cols)[self.usecols]:
                    if not isinstance(dim, list):
                        dim = [dim]
                    yield df1, dim
            else:
                while input_handle.shape[0] > 0:
                    df = input_handle.iloc[:self.chunksize]
                    for dim, df1 in df.groupby(self.dim_cols)[self.usecols]:
                        if not isinstance(dim, list):
                            dim = [dim]
                        yield df1, dim
                    input_handle = input_handle.iloc[self.chunksize:]
        else:
            yield from _input_parser(input_handle, self.Formats, self.sep,
                                     self.usecols, self.dim_cols,
                                     self.chunksize)

    def tocz(self, Input=None, usecols=[4, 5], dim_cols=[0],
             sep='\t', chunksize=5000, header=None, skiprows=0):
        """
        Pack dataframe, stdin or a file path into .cz file with or without reference
        coordinates file.

        Parameters
        ----------
        Input : input
            list, tuple, np.ndarray or dataframe, stdin (Input is None, "stdin" or "-"),
            or a file path (need to specify sep,header and skiprows).
        usecols : list
            usecols is the index of columns to be packed into .cz Columns.
        dim_cols : list
            index of columns to be set as Dimensions of Writer, such as chrom
            columns.
        sep : str
            defauls is "\t", used as separator for the input file path.
        chunksize : int
            Chunk Input dataframe of file path.
        header : bool
            whether the input file path contain header.
        skiprows : int
            number of rows to skip for input file path.
        pr: int or str
            If reference was given, use the `pr` column as coordinates, which
            would be used to match the `p` column in the Input file or df.
            If Input is a file path, then `pr` should be int, if Input is a dataframe with
            custom columns names, then `pr` should be string matching the column
            name in dataframe.
        p: int or str
            Similar as `pr`, `p` is the column in Input to match the coordinate in
            the reference file. When Input is a file path or stdin, p should be int,
            but if Input is a dataframe, p should be the name appear in the dataframe
            header columns.

        Returns
        -------

        """
        self.sep = sep
        if not isinstance(Input,pd.DataFrame):
            if isinstance(usecols, int):
                self.usecols = [int(usecols)]
            elif isinstance(usecols, str):
                self.usecols = [int(i) for i in usecols.split(',')]
            else:
                self.usecols = [int(i) for i in usecols]
            if isinstance(dim_cols, int):
                self.dim_cols = [int(dim_cols)]
            elif isinstance(dim_cols, str):
                self.dim_cols = [int(i) for i in dim_cols.split(',')]
            else:
                self.dim_cols = [int(i) for i in dim_cols]
        else:
            self.usecols = usecols
            self.dim_cols = dim_cols

        assert len(self.usecols) == len(self.Formats)
        assert len(self.dim_cols) == len(self.Dimensions)
        self.chunksize = chunksize
        self.header = header
        self.skiprows = skiprows
        # if self.verbose > 0:
        #     print("Input: ", type(Input), Input)

        if Input is None or Input == 'stdin':
            data_generator = self._parse_input_no_ref(sys.stdin)
        elif isinstance(Input, str):
            input_path = os.path.abspath(os.path.expanduser(Input))
            if not os.path.exists(input_path):
                raise ValueError("Unknown format for Input")
            if os.path.exists(input_path + '.tbi'):
                print("Please use bed2cz command to convert input to .cz.")
                return
            else:
                data_generator = self._parse_input_no_ref(input_path)
        elif isinstance(Input, (list, tuple, np.ndarray)):
            Input = pd.DataFrame(Input)
            data_generator = self._parse_input_no_ref(Input)
        elif isinstance(Input, pd.DataFrame):
            data_generator = self._parse_input_no_ref(Input)
        else:
            raise ValueError(f"Unknow input type")

        # if self.verbose > 0:
        #     print(self.usecols, self.dim_cols)
        #     print(type(self.usecols), type(self.dim_cols))
        for df, dim in data_generator:
            data = df.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                            axis=1).sum()
            self.write_chunk(data, dim)
        self.close()

    @staticmethod
    def create_new_dim(basename):
        if basename.endswith('.cz'):
            return basename[:-3]
        return basename

    def catcz(self, Input=None, dim_order=None, add_dim=False,
              title="filename"):
        """
		Cat multiple .cz files into one .cz file.

		Parameters
		----------
		Input : str or list
			Either a str (including *, as input for glob, should be inside the\
			double quotation marks if using fire) or a list.
		dim_order : None, list or path
			If dim_order=None, Input will be sorted using python sorted.
			If dim_order is a list, tuple or array of basename.rstrip(.cz), sorted as dim_order.
			If dim_order is a file path (for example, chrom size path to dim_order chroms
			or only use selected chroms) will be sorted as
			the 1st column of the input file path (without header, tab-separated).
			default is None
        add_dim: bool or function
            whether to add .cz file names as an new dimension to the merged
            .cz file. For example, we have multiple .cz files for many cells, in
            each .cz file, the dimensions are ['chrom'], after merged, we would
            like to add file name of .cz as a new dimension ['cell_id']. In this case,
            the dimensions in the merged header would be ["chrom","cell_id"], and
            in each chunk, in addition to the previous dim ['chr1'] or ['chr22'].., a
            new dim would also be append to the previous dim, like ['chr1','cell_1'],
            ['chr22','cell_100'].
            However, if add_dim is a function, the input to this function is the .cz
            file basename, the returned value from this funcion would be used
            as new dim and added into the chunk_dims. The default function to
            convert filename to dim name is self.create_new_dim.
        title: str
            if add_dim is True or a python function, title would be append to
            the header['Dimensions'] of the merged .cz file's header. If the title of
            new dimension had already given in Writer Dimensions,
            title can be None, otherwise, title should be provided.

		Returns
		-------

		"""
        if isinstance(Input, str) and '*' in Input:
            Input = glob.glob(Input)
        if not isinstance(Input, list):
            raise ValueError("Input should be either a list of a string including *.")
        if dim_order is None:
            Input = sorted(Input)
        else:
            # creating a dict, keys are file base name, value are full path
            D = {os.path.basename(inp)[:-3]: inp for inp in Input}
            if self.verbose > 0:
                print(D)
            if isinstance(dim_order, str):
                # read each file from path containing file basename
                dim_order = pd.read_csv(os.path.abspath(os.path.expanduser(dim_order)),
                                    sep='\t', header=None, usecols=[0])[0].tolist()
            if isinstance(dim_order, (list, tuple, np.ndarray)):  # dim_order is a list
                # Input=[str(i)+'.cz' for i in dim_order]
                Input = [D[str(i)] for i in dim_order]
            else:
                raise ValueError("input of dim_order is not corrected !")
        if self.verbose > 0:
            print(Input)
        self.new_dim_creator = None
        if add_dim!=False: # add filename as another dimension.
            if add_dim==True:
                self.new_dim_creator=self.create_new_dim
            elif callable(add_dim): #is a function
                self.new_dim_creator = add_dim
            else:
                raise ValueError("add_dim should either be True,False or a function")

        for file_path in Input:
            reader = Reader(file_path)
            # data_size = reader.header['total_size'] - reader.header['header_size']
            chunks = reader.get_chunks()
            (start_offset, chunk_size, dims, data_len,
             end_offset, nblocks, b1str_virtual_offsets) = next(chunks)
            # check whether new dim has already been added to header
            if not self.new_dim_creator is None:
                new_dim = [self.new_dim_creator(os.path.basename(file_path))]
            else:
                new_dim = []
            if len(dims + tuple(new_dim)) > len(self.Dimensions):
                # new_dim title should be also added to header Dimensions
                self.Dimensions = self.Dimensions + [title]
                self._handle.seek(self._n_dim_offset)
                # when a new dim is added, go back here and rewrite the n_dim (1byte)
                self._handle.write(struct.pack("<B", len(self.Dimensions)))  # 1byte
                self._handle.seek(self._header_size)
                # write new dim to the end of header
                dname_len = len(title)
                self._handle.write(struct.pack("<B", dname_len))  # length of each name, 1 byte
                self._handle.write(struct.pack(f"<{dname_len}s", bytes(title, 'utf-8')))
                self._header_size = self._handle.tell()
            while True:  # process each file (chunk)
                reader._handle.seek(start_offset)  # chunk start offset
                # write the enrire chunk
                new_chunk_start_offset = self._handle.tell()
                delta_offset = new_chunk_start_offset - start_offset
                # self._handle.write(reader._handle.read(end_offset - start_offset))
                self._handle.write(reader._handle.read(chunk_size + 16))
                # modify the chunk_black_1st_record_virtual_offsets
                new_block_1st_record_vof = []
                for offset in b1str_virtual_offsets:
                    block_offset, within_block_offset = split_virtual_offset(offset)
                    new_block_offset = delta_offset + block_offset
                    new_1st_rd_vof = make_virtual_offset(new_block_offset,
                                                         within_block_offset)
                    new_block_1st_record_vof.append(new_1st_rd_vof)
                # write new_block_1st_record_vof
                for block_offset in new_block_1st_record_vof:
                    self._handle.write(struct.pack("<Q", block_offset))
                # rewrite the new dim onto the tail of chunk
                for dim in dims + tuple(new_dim):
                    dim_len = len(dim)
                    self._handle.write(struct.pack("<B", dim_len))  # length of each dim, 1 byte
                    self._handle.write(struct.pack(f"<{dim_len}s", bytes(dim, 'utf-8')))
                    # print(dim_len,dim)
                # read next chunk
                try:
                    (start_offset, chunk_size, dims, data_len,
                     end_offset, nblocks, b1str_virtual_offsets) = chunks.__next__()
                except:
                    break
            reader.close()
        self._write_total_size_eof_close()

    def flush(self):
        """Flush data explicitally."""
        while len(self._buffer) >= _BLOCK_MAX_LEN:
            # I don't think this is going to happen
            self._write_block(self._buffer[:_BLOCK_MAX_LEN])
            self._buffer = self._buffer[_BLOCK_MAX_LEN:]
        self._write_block(self._buffer)
        self._chunk_finished()
        self._buffer = b""
        self._handle.flush()

    def _write_total_size_eof_close(self):
        cur_offset = self._handle.tell()
        self._handle.seek(self._magic_size + 4)  # magic and version
        self._handle.write(struct.pack("<Q", cur_offset))  # real offset.
        self._handle.seek(cur_offset)
        self._handle.write(_bcz_eof)
        self._handle.flush()
        self._handle.close()

    def close(self):
        """Flush data, write 28 bytes BGZF EOF marker, and close BGZF file.

		samtools will look for a magic EOF marker, just a 28 byte empty BGZF
		block, and if it is missing warns the BAM file may be truncated. In
		addition to samtools writing this block, so too does bgzip - so this
		implementation does too.
		"""
        if self._buffer:
            self.flush()
        else:
            self._chunk_finished()
        self._write_total_size_eof_close()

    def tell(self):
        """Return a BGZF 64-bit virtual offset."""
        return make_virtual_offset(self._handle.tell(), len(self._buffer))

    def seekable(self):
        return False

    def isatty(self):
        """Return True if connected to a TTY device."""
        return False

    def fileno(self):
        """Return integer file descriptor."""
        return self._handle.fileno()

    def __enter__(self):
        """Open a file operable with WITH statement."""
        return self

    def __exit__(self, type, value, traceback):
        """Close a file with WITH statement."""
        self.close()
# ==========================================================
def test_difference(file1, file2, usecols1, usecols2, sep='\t', header1=None,
                    header2=None, intersect=False, index_cols=[0, 1]):
    df1 = pd.read_csv(file1, usecols=usecols1, sep=sep, header=header1)
    df2 = pd.read_csv(file2, usecols=usecols2, sep=sep, header=header2)
    df1.columns = list(range(df1.shape[1]))
    df2.columns = list(range(df2.shape[1]))
    if intersect:
        df1.set_index(index_cols, inplace=True)
        df2.set_index(index_cols, inplace=True)
        common_index = list(set(df1.index.tolist()) & set(df2.index.tolist()))
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]
    else:
        df1.sort_values(usecols1, inplace=True)
        df2.sort_values(usecols2, inplace=True)
    for i in range(min(len(usecols1), len(usecols2))):
        print(i)
        print(df1.iloc[:, i].tolist() == df2.iloc[:,i].tolist())

# ==========================================================

if __name__=="__main__":
    import fire
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire()