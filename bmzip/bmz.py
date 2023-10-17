#!/usr/bin/env python
import glob
import os, sys
import struct
import zlib
from builtins import open as _open
import numpy as np
import pandas as pd
import gzip

_bmz_magic = b'BMZIP'
_block_magic = b"MB"
_chunk_magic = b"MC"
_BLOCK_MAX_LEN = 65535
_bmz_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BM\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
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
def get_dtfuncs(formats):
    return [dtype_func[t[-1]] for t in formats]


# ==========================================================
def make_virtual_offset(block_start_offset, within_block_offset):
    """Compute a BGZF virtual offset from block start and within block offsets.

	The BAM indexing scheme records read positions using a 64 bit
	'virtual offset', comprising in C terms:

	block_start_offset << 16 | within_block_offset

	Here block_start_offset is the file offset of the BGZF block
	start (unsigned integer using up to 64-16 = 48 bits), and
	within_block_offset within the (decompressed) block (unsigned
	16 bit integer).
	"""
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
    """Divides a 64-bit BGZF virtual offset into block start & within block offsets.

	>>> (100000, 0) == split_virtual_offset(_BLOCK_MAX_LEN00000)
	True
	>>> (100000, 10) == split_virtual_offset(_BLOCK_MAX_LEN00010)
	True

	"""
    start = virtual_offset >> 16
    return start, virtual_offset ^ (start << 16)


# ==========================================================
def SummaryBmzBlocks(handle):
    """Low level debugging function to inspect BGZF blocks.

	Expects a BGZF compressed file opened in binary read mode using
	the builtin open function. Do not use a handle from this bgzf
	module or the gzip module's open function which will decompress
	the file.

	Returns the block start offset (see virtual offsets), the block
	length (add these for the start of the next block), and the
	decompressed length of the blocks contents (limited to _BLOCK_MAX_LEN in
	BGZF), as an iterator - one tuple per BGZF block.

	>>> from builtins import open
	>>> handle = open("SamBam/ex1.bam", "rb")
	>>> for values in BmzBlocks(handle):
	...     print("Raw start %i, raw length %i; data start %i, data length %i" % values)
	Raw start 0, raw length 18239; data start 0, data length _BLOCK_MAX_LEN
	Raw start 18239, raw length 18223; data start _BLOCK_MAX_LEN, data length _BLOCK_MAX_LEN
	Raw start 36462, raw length 18017; data start 131072, data length _BLOCK_MAX_LEN
	Raw start 54479, raw length 17342; data start 196608, data length _BLOCK_MAX_LEN
	Raw start 71821, raw length 17715; data start 262144, data length _BLOCK_MAX_LEN
	Raw start 89536, raw length 17728; data start 327680, data length _BLOCK_MAX_LEN
	Raw start 107264, raw length 17292; data start 393216, data length 63398
	Raw start 124556, raw length 28; data start 456614, data length 0
	>>> handle.close()

	Indirectly we can tell this file came from an old version of
	samtools because all the blocks (except the final one and the
	dummy empty EOF marker block) are _BLOCK_MAX_LEN bytes.  Later versions
	avoid splitting a read between two blocks, and give the header
	its own block (useful to speed up replacing the header). You
	can see this in ex1_refresh.bam created using samtools 0.1.18:

	samtools view -b ex1.bam > ex1_refresh.bam

	>>> handle = open("SamBam/ex1_refresh.bam", "rb")
	>>> for values in BmzBlocks(handle):
	...     print("Raw start %i, raw length %i; data start %i, data length %i" % values)
	Raw start 0, raw length 53; data start 0, data length 38
	Raw start 53, raw length 18195; data start 38, data length 65434
	Raw start 18248, raw length 18190; data start 65472, data length 65409
	Raw start 36438, raw length 18004; data start 130881, data length 65483
	Raw start 54442, raw length 17353; data start 196364, data length 65519
	Raw start 71795, raw length 17708; data start 261883, data length 65411
	Raw start 89503, raw length 17709; data start 327294, data length 65466
	Raw start 107212, raw length 17390; data start 392760, data length 63854
	Raw start 124602, raw length 28; data start 456614, data length 0
	>>> handle.close()

	The above example has no embedded SAM header (thus the first block
	is very small at just 38 bytes of decompressed data), while the next
	example does (a larger block of 103 bytes). Notice that the rest of
	the blocks show the same sizes (they contain the same read data):

	>>> handle = open("SamBam/ex1_header.bam", "rb")
	>>> for values in BmzBlocks(handle):
	...     print("Raw start %i, raw length %i; data start %i, data length %i" % values)
	Raw start 0, raw length 104; data start 0, data length 103
	Raw start 104, raw length 18195; data start 103, data length 65434
	Raw start 18299, raw length 18190; data start 65537, data length 65409
	Raw start 36489, raw length 18004; data start 130946, data length 65483
	Raw start 54493, raw length 17353; data start 196429, data length 65519
	Raw start 71846, raw length 17708; data start 261948, data length 65411
	Raw start 89554, raw length 17709; data start 327359, data length 65466
	Raw start 107263, raw length 17390; data start 392825, data length 63854
	Raw start 124653, raw length 28; data start 456679, data length 0
	>>> handle.close()

	"""
    if isinstance(handle, Reader):
        raise TypeError("Function BmzBlocks expects a binary handle")
    data_start = 0
    while True:
        start_offset = handle.tell()
        try:
            block_length, data_len = _load_bmz_block(handle)
        except StopIteration:
            break
        yield start_offset, block_length, data_start, data_len
        data_start += data_len


# ==========================================================
def _load_bmz_block(handle, decompress=False):
    """Load the next BGZF block of compressed data (PRIVATE).

	Returns a tuple (block size and data), or at end of file
	will raise StopIteration.
	"""
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
    if infile.endswith('.gz'):
        f = gzip.open(infile, 'rb')
    else:
        f = open(infile, 'r')
    return f


# ==========================================================
class Reader:
    def __init__(self, Input=None, mode="rb", fileobj=None, max_cache=100):
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
        if magic != _bmz_magic:
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

    def print_header(self):
        for k in self.header:
            print(k," : ",self.header[k])

    def _load_chunk(self, start_offset=None, view=True):
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
        if view:  # no need to load _chunk_block_1st_record_virtual_offsets
            self._handle.seek(self._chunk_nblocks * 8, 1)
        # _chunk_tail_offset = end position of this chunk
        # = start position of next chunk.
        else:
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
        r = self._load_chunk(self.header['header_size'], view=False)
        while r:
            yield [self._chunk_start_offset, self._chunk_size,
                   self._chunk_dims, self._chunk_data_len, self._chunk_end_offset,
                   self._chunk_nblocks, self._chunk_block_1st_record_virtual_offsets
                   ]
            r = self._load_chunk(view=False)

    def summary_chunks(self,print=True):
        r = self._load_chunk(self.header['header_size'], view=True)
        header = ['chunk_start_offset', 'chunk_size', 'chunk_dims', 'chunk_data_len',
                  'chunk_tail_offset', 'chunk_nblocks']
        if print:
            sys.stdout.write('\t'.join(header) + '\n')
        else:
            R = []
        while r:
            self._handle.seek(self._chunk_start_offset + 10)
            chunk_info = [self._chunk_start_offset, self._chunk_size,
                          self._chunk_dims, self._chunk_data_len, self._chunk_end_offset,
                          self._chunk_nblocks]
            try:
                if print:
                    sys.stdout.write('\t'.join([str(v) for v in chunk_info]) + '\n')
                else:
                    R.append(chunk_info)
            except:
                sys.stdout.close()
            r = self._load_chunk(view=True)
        if print:
            sys.stdout.close()
        else:
            df = pd.DataFrame(R, columns=header)
            return df

    def summary_blocks(self, print=True):
        r = self._load_chunk(self.header['header_size'], view=True)
        header = ['chunk_dims'] + ['block_start_offset','block_size',
                                   'block_data_start']
        if print:
            sys.stdout.write('\t'.join(header) + '\n')
        else:
            R=[]
        while r:
            self._handle.seek(self._chunk_start_offset + 10)
            chunk_info = [self._chunk_dims]
            for block in SummaryBmzBlocks(self._handle):
                block = chunk_info + list(block)[:-1]
                if print:
                    try:
                        sys.stdout.write('\t'.join([str(v) for v in block]) + '\n')
                    except:
                        sys.stdout.close()
                else:
                    R.append(block)
            r = self._load_chunk(view=True)
        if print:
            sys.stdout.close()
        else:
            df=pd.DataFrame(R,columns=header)
            return df

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
            block_size, self._buffer = _load_bmz_block(self._handle, True)
        except StopIteration:  # EOF
            block_size = 0
            self._buffer = b""
        self._within_block_offset = 0
        self._block_raw_length = block_size

    def byte2str(self,values):
        return [str(v,'utf-8') if f[-1] in ['s','c'] else str(v)
                for v,f in zip(values,self.header['Formats'])]

    def _print_cache(self):
        self._cached_data += self._buffer
        end_index = len(self._cached_data) - (len(self._cached_data) % self._unit_size)
        for result in struct.iter_unpack(f"<{self.fmts}", self._cached_data[:end_index]):
            line = '\t'.join(self.byte2str(result))
            try:
                sys.stdout.write(self.chunk_dim + line + '\n')
            except:
                sys.stdout.close()
        self._cached_data = self._cached_data[end_index:]

    def view(self, show_dim=None, header=True, dim=None):
        """
		View .bmz file.

		Parameters
		----------
		show_dim : str
			index of dims given to writer.write_chunk, separated by comma,
			default is None, dims[show_dim] will be shown in each row (such as sampleID
			and chrom)
		header : bool
			whether to print the header.
        dim: None, bool, list or file path
            None (default): use the default order in .mz file;

            bool (True): sort the dims and used as dim

            list: use this dim (dims) as order and print records.

            path: use the first len(dim) columns as dim order, there should be
                no header in file path and use \t as separator.
		Returns
		-------

		"""
        if type(show_dim) == int:
            show_dim = [show_dim]
        elif type(show_dim)==str:
            show_dim = [int(i) for i in show_dim.split(',')]

        chunk_info = self.summary_chunks(print=False)
        chunk_info.set_index('chunk_dims', inplace=True) #chunk_dims is a tuple.
        if not dim is None:
            # equal to query chromosome if self.header['dimensions'][0]==chrom
            if isinstance(dim, str):
                order_path=os.path.abspath(os.path.expanduser(dim))
                if os.path.exists(order_path):
                    dim = pd.read_csv(order_path,sep='\t', header=None,
                                    usecols=show_dim)[show_dim].apply(lambda x:tuple(x.tolist()),
                                                            axis=1).tolist()
                else:
                    dim=dim.split(',')
            if type(dim[0]) == str: # each element of dim should be a tuple
                dim = [tuple([o]) for o in dim]
            if not isinstance(dim, (list, tuple, np.ndarray)):  # dim is a list
                raise ValueError("input of order is not corrected !")
            chunk_info=chunk_info.loc[dim]

        if not show_dim is None:
            header_dim = "\t".join([self.header['Dimensions'][t] for t in show_dim]) + '\t'
        else:
            header_dim = ''
        if header: #show header
            line = "\t".join(self.header['Columns'])
            sys.stdout.write(header_dim + line + '\n')

        for start_pos in chunk_info.chunk_start_offset.tolist():
            r = self._load_chunk(start_pos, view=True) #self.header['header_size']
            # header_size is the start offset of 1st chunk
            self._cached_data = b''
            # here, we have already read the chunk header and tail but
            # not the blocks, we need to jump to body start pos and read each block.
            if not show_dim is None:
                self.chunk_dim = "\t".join([self._chunk_dims[t] for t in show_dim]) + '\t'
            else:
                self.chunk_dim = ''
            self._load_block(start_offset=self._chunk_start_offset + 10)  #
            while self._block_raw_length > 0:
                # deal with such case: unit_size is 10, but data(_buffer) size is 18,
                self._print_cache()
                self._load_block()
            # r = self._load_chunk(view=True)
        sys.stdout.close()

    def _read_1record(self):
        tmp=self._buffer1[:self._unit_size]
        self._buffer1 = self._buffer1[self._unit_size: ]
        return struct.unpack(f"<{self.fmts}", tmp)

    def _query(self, regions,s,e):
        chunk_info = self.summary_chunks(print=False)
        dim2pos=chunk_info.set_index('chunk_dims').chunk_start_offset.to_dict()
        for dim,start,end in regions:
            r = self._load_chunk(dim2pos[dim], view=False)
            # find the closest block to a given start position
            block_1st_records=[
                self._seek_and_read_1record(offset)
                for offset in self._chunk_block_1st_record_virtual_offsets
            ]
            block_1st_starts=[r[s] for r in block_1st_records]
            start_block_index=0
            while True:
                if start_block_index+1 >= self._chunk_nblocks:
                    break
                if block_1st_starts[start_block_index+1] > start:
                    break
                start_block_index+=1
            if s == e:
                end_block_index = start_block_index
            else:
                block_1st_ends = [r[e] for r in block_1st_records]
                end_block_index = 0
                while True:
                    if end_block_index +1 >= self._chunk_nblocks:
                        break
                    if block_1st_ends[end_block_index+1] > end:
                        break
                    end_block_index += 1
            for block_index in range(start_block_index,end_block_index+1):
                # print(block_index)
                virtual_offset=self._chunk_block_1st_record_virtual_offsets[block_index]
                self.seek(virtual_offset) #seek to the target block, self._buffer
                self._buffer1=self._buffer[self._within_block_offset:]
                record=self._read_1record()
                while record[s] < start:
                    if len(self._buffer1) < self._unit_size:
                        break
                    record = self._read_1record()
                while record[e] <= end:
                    # print(list(dim) + list(record))
                    yield dim,record
                    if len(self._buffer1) < self._unit_size:
                        break
                    record = self._read_1record()

    def query(self, dim=None,start=None,end=None,regions=None,
               query_col=[0],print=True):
        """
        query .mz file by dim, start and end, if regions provided, dim, start and
        end should be None, regions should be a list, each element of regions
        is a list, for example, regions=[[('cell1','chr1'),1,10],
        [('cell10','chr22'),100,200]],and so on.

        Parameters
        ----------
        dim : str
            dimension, such as chr1 (if Dimension 'chrom' is inlcuded in
            header['Dimension']), or sample1 (if something like 'sampleID' is
            included in header['Dimension'] and chunk contains such dimension)
        start : int
            start position, if None, the entire dim would be returned.
        end : int
            end position
        regions : list or file path
            regions=[
                [dim,start,end], [dim,start,end]
            ],
            dim is a tuple, for example, dim=('chr1');

            if regions is a file path (separated by tab and no header),
            the first len(header['Dimensions']) columns would be
            used as dim, and next two columns would be used as start and end.
        query_col: list
            index of columns (header['Columns']) to be queried,for example,
            if header['Columns']=['pos','mv','cov'], then pos is what we want to
            query, so query_col should be [0], but if
            header['Columns']=['start','end','peak'], then start and end are what
            we would like to query, query_col should be [0,1]

        Returns
        -------

        """
        # dim,start,end=('chr12',),3109883,3110041
        if dim is None and regions is None:
            raise ValueError("Please provide either dim,start,end or regions")
        if (not dim is None) and (not regions is None):
            raise ValueError("Please query by dim,start,end or by regions, "
                             "not both !")
        if not dim is None:
            if type(dim)==str:
                dim=tuple([dim])
            regions=[[dim,start,end]]
        else:
            if isinstance(regions,str):
                region_path=os.path.abspath(os.path.expanduser(regions))
                if os.path.exists(region_path):
                    n_dim=self.header['Dimensions']
                    usecols=list(range(n_dim+2))
                    df=pd.read_csv(region_path,sep='\t',header=None,
                                   usecols=usecols)
                    regions=df.apply(lambda x:[
                        tuple(x[:n_dim]),
                        x[n_dim],
                        x[n_dim+1]
                    ],axis=1).tolist()
                else:
                    raise ValueError(f"path {region_path} not existed.")
            else: #check format of regions
                if type(regions[0][0]) != tuple:
                    raise ValueError("The first element of first region is not a tuple:"
                                     f"{type(regions[0][0])}")

        if len(query_col)==1:
            s = e = query_col[0]
        elif len(query_col)==2:
            s,e=query_col
        else:
            raise ValueError("length of query_col can not be greater then 2.")

        header=self.header['Dimensions'] + self.header['Columns']
        if print:
            sys.stdout.write('\t'.join(header)+'\n')
            for dims,row in self._query(regions,s,e):
                sys.stdout.write('\t'.join([str(d) for d in dims]+self.byte2str(row))+'\n')
            sys.stdout.close()
        else:
            for dims,row in self._query(regions,s,e):
                yield list(dims)+list(row)

    def _seek_and_read_1record(self,virtual_offset):
        self.seek(virtual_offset)
        return struct.unpack(f"<{self.fmts}", self.read(self._unit_size))

    def build_index(self, cols=[0]):
        """
        create index for the .mz file passed to the Reader.

        Parameters
        ----------
        name_cols : int, list, or str
            name_cols (default is [0]) should be the index in the names list (
            see self.header['names'], or bmzip Reader -I test.mz print_header)
            if a int given, build index for only one name,
            if a list give, build index for multiple names,
            if a str give, the string will be convert to list (split by comma)

        Returns
        -------

        """
        if type(cols) == int:
            cols=[cols]
        elif type(cols) ==str:
            cols = [int(i) for i in cols.split(',')]
        else:
            cols=[int(i) for i in cols]

        block_info=self.summary_blocks(print=False)
        n_chunk_dims=len(block_info.chunk_dims.iloc[0])
        Dimension=[]
        for i,dim in zip(range(n_chunk_dims),self.header['Dimensions']):
            # print(i,dim)
            Dimension.append(dim)
            block_info[dim]=block_info.chunk_dims.apply(lambda x:x[i])
        block_info['within_block_offset'] = (
            np.ceil(block_info.block_data_start / self._unit_size) * self._unit_size
            - block_info.block_data_start
        ).map(int)
        block_info['1st_record_virtual_offset'] = (
            block_info.apply(lambda x:
                make_virtual_offset(x.block_start_offset, x.within_block_offset),
                             axis=1)
        )

        block_info['1st_record']=block_info['1st_record_virtual_offset'].apply(
            self._seek_and_read_1record)

        columns = []
        for col in cols:
            Col = self.header['Columns'][col]
            block_info[Col] = block_info['1st_record'].apply(lambda x:x[col])
            columns.append(Col)
        formats=[self.header['Formats'][i] for i in cols]
        self.bmi=self.Input + '.bmi'
        Columns = columns+['1st_record_virtual_offset']
        self.bmi_writer=Writer(Output=self.bmi, Formats=formats+['Q'],
                               Columns=Columns,Dimension=Dimension)
        self.bmi_writer.pack(Input=block_info.loc[:,Dimension+Columns],
                             usecols=Columns,dim_cols=Dimension,
                             chunksize=None)

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

    def read(self, size=-1):
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
                if not data:
                    raise ValueError("Must be at least 1 byte")
                result += data
                break
            else:
                data = self._buffer[self._within_block_offset:]
                size -= len(data)
                self._load_block()  # will reset offsets
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

    def __next__(self):
        """Return the next line."""
        line = self.readline()
        if not line:
            raise StopIteration
        return line

    def __iter__(self):
        """Iterate over the lines in the BGZF file."""
        return self

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
class Writer:
    def __init__(self, Output=None, mode="wb", Formats=['H', 'H'],
                 Columns=['mc', 'cov'], Dimension=['chrom'], fileobj=None,
                 message='',level=6, verbose=0):
        """
        bmzip .mz writer.

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
        Dimension : list
            Dimension to be included for each chunk, for example, if you would like
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
                    raise ValueError(f"Must use write or append mode, not {mode!r}")
                handle = _open(Output, mode)
            else:  # write to stdout buffer
                handle = sys.stdout.buffer
        self._handle = handle
        self._buffer = b""
        self._chunk_start_offset = None
        self._chunk_dims = None
        if type(Formats) == str and ',' not in Formats:
            self.Formats = [Formats]
        elif type(Formats) == str:
            self.Formats = Formats.split(',')
        else:
            self.Formats = list(Formats)
        self.level = level
        if type(Columns) == str:
            self.Columns = Columns.split(',')
        else:
            self.Columns = list(Columns)
        if type(Dimension) == str:
            self.Dimension = Dimension.split(',')
        else:
            self.Dimension = list(Dimension)
        self._magic_size = len(_bmz_magic)
        self.verbose = verbose
        self.message=message
        if self.verbose > 0:
            print(self.Formats, self.Columns, self.Dimension)
            print(type(self.Formats), type(self.Columns), type(self.Dimension))
        self.write_header()

    def write_header(self):
        f = self._handle
        f.write(struct.pack(f"<{self._magic_size}s", _bmz_magic))  # 5 bytes
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
        f.write(struct.pack("<B", len(self.Dimension)))  # 1byte
        for dim in self.Dimension:
            dname_len = len(dim)
            f.write(struct.pack("<B", dname_len))  # length of each name, 1 byte
            f.write(struct.pack(f"<{dname_len}s", bytes(dim, 'utf-8')))
        # when multiple .mz are cat into one .mz and a new dimension is added,
        # such as sampleID, seek to this position (_header_size) and write another
        # two element: new_dname_len (B) and new_dim, then go to
        # _n_dim_offset,rewrite the n_dim (n_dim = n_dim + 1) and
        # self.Dimension.append(new_dim)
        self._header_size=f.tell()
        self.fmts = ''.join(list(self.Formats))
        self._unit_size = struct.calcsize(self.fmts)

    def _write_block(self, block):
        if len(block) > _BLOCK_MAX_LEN:  # 65536 = 1 << 16 = 2**16
            raise ValueError(f"{len(block)} Block length > {_BLOCK_MAX_LEN}")
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
        if len(compressed) > _BLOCK_MAX_LEN:
            raise RuntimeError(
                "TODO - Didn't compress enough, try less data in this block"
            )
        bsize = struct.pack("<H", len(compressed) + 6)
        # block size: magic (2 btyes) + block_size (2bytes) + compressed data +
        # block_data_len (2 bytes)
        uncompressed_length = struct.pack("<H", len(block))  # 2 bytes
        data = _block_magic + bsize + compressed + uncompressed_length
        block_start_offset=self._handle.tell()
        within_block_offset = int(
            np.ceil(self._chunk_data_len / self._unit_size) * self._unit_size
            - self._chunk_data_len)
        virtual_offset=(block_start_offset << 16) | within_block_offset
        self._block_1st_record_virtual_offsets.append(virtual_offset)
        # _block_offsets are real position on disk, not a virtual offset
        self._handle.write(data)
        self._chunk_data_len += len(block)

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
            as Dimension given to Writer.

        Returns
        -------

        """
        assert len(dims)==len(self.Dimension)
        if isinstance(data, str):
            data = data.encode("latin-1")
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
            if self.verbose > 0:
                print("Writing chunk with dims: ", self._chunk_dims)

        if len(self._buffer) + len(data) < _BLOCK_MAX_LEN:
            self._buffer += data
        else:
            self._buffer += data
            while len(self._buffer) >= _BLOCK_MAX_LEN:
                self._write_block(self._buffer[:_BLOCK_MAX_LEN])
                self._buffer = self._buffer[_BLOCK_MAX_LEN:]

    def parse_input_pd(self, input_handle):
        if isinstance(input_handle, pd.DataFrame):
            # usecols and dim_cols should be in the columns of this dataframe.
            if self.chunksize is None:
                for dim, df1 in input_handle.groupby(self.dim_cols)[self.usecols]:
                    if type(dim) != list:
                        dim = [dim]
                    yield df1.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                                    axis=1).sum(), dim
            else:
                while input_handle.shape[0] > 0:
                    df = input_handle.iloc[:self.chunksize]
                    for dim, df1 in df.groupby(self.dim_cols)[self.usecols]:
                        if type(dim) != list:
                            dim = [dim]
                        yield df1.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                                        axis=1).sum(), dim
                    input_handle = input_handle.iloc[self.chunksize:]
        else:  # stdin or read from file.
            for df in pd.read_csv(input_handle, sep=self.sep,
                                  usecols=self.dim_cols + self.usecols,
                                  chunksize=self.chunksize, header=self.header,
                                  skiprows=self.skiprows):
                for dim, df1 in df.groupby(self.dim_cols)[self.usecols]:
                    if type(dim) != list:
                        dim = [dim]
                    yield df1.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                                    axis=1).sum(), dim

    def parse_input(self, input_handle):
        if isinstance(input_handle, pd.DataFrame):
            # usecols and dim_cols should be in the columns of this dataframe.
            if self.chunksize is None:
                for dim, df1 in input_handle.groupby(self.dim_cols)[self.usecols]:
                    if type(dim) != list:
                        dim = [dim]
                    yield df1.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                                    axis=1).sum(), dim
            else:
                while input_handle.shape[0] > 0:
                    df = input_handle.iloc[:self.chunksize]
                    for dim, df1 in df.groupby(self.dim_cols)[self.usecols]:
                        if type(dim) != list:
                            dim = [dim]
                        yield df1.apply(lambda x: struct.pack(f"<{self.fmts}", *x.tolist()),
                                        axis=1).sum(), dim
                    input_handle = input_handle.iloc[self.chunksize:]
        else:
            if hasattr(input_handle, 'readline'):  # stdin or read from file.
                f = input_handle
            else:  # input_handle is a file path
                f = open1(input_handle)
            data, i, pre_dims = b'', 0, None
            dtfuncs = get_dtfuncs(self.Formats)
            line = f.readline()
            while line:
                if isinstance(line, bytes):
                    line = line.decode('utf-8')
                values = line.rstrip('\n').split(self.sep)
                dims = [values[i] for i in self.dim_cols]
                if dims != pre_dims:  # a new dim (chrom), for example, chr1 -> chr2
                    if len(data) > 0:  # write rest data of chr1
                        yield data, pre_dims
                        data, i = b'', 0
                    pre_dims = dims
                if i >= self.chunksize:  # dims are the same, but reach chunksize
                    yield data, pre_dims
                    data, i = b'', 0
                values_to_pack = [values[i] for i in self.usecols]
                data += struct.pack(f"<{self.fmts}",
                                    *[func(v) for v, func in zip(values_to_pack, dtfuncs)]
                                    )
                line = f.readline()
                i += 1
            f.close()

    def pack(self, Input=None, usecols=[1, 4, 5], dim_cols=[0],
             sep='\t', chunksize=5000, header=None, skiprows=0):
        self.sep = sep
        if not isinstance(Input,pd.DataFrame):
            if type(usecols) == int:
                self.usecols = [int(usecols)]
            elif type(usecols) == str:
                self.usecols = [int(i) for i in usecols.split(',')]
            else:
                self.usecols = [int(i) for i in usecols]
            if type(dim_cols) == int:
                self.dim_cols = [int(dim_cols)]
            elif type(dim_cols) == str:
                self.dim_cols = [int(i) for i in dim_cols.split(',')]
            else:
                self.dim_cols = [int(i) for i in dim_cols]
        else:
            self.usecols=usecols
            self.dim_cols=dim_cols

        assert len(self.usecols) == len(self.Formats)
        assert len(self.dim_cols) == len(self.Dimension)
        self.chunksize = chunksize
        self.header = header
        self.skiprows = skiprows
        if isinstance(Input, (list, tuple, np.ndarray)):
            Input = pd.DataFrame(Input)
        if isinstance(Input, str):
            input_path = os.path.abspath(os.path.expanduser(Input))
            data_generator = self.parse_input(input_path)
        elif isinstance(Input,pd.DataFrame):  # Input is a dataframe
            data_generator = self.parse_input(Input)
        elif Input is None or Input == 'stdin' or Input == '-':
            data_generator = self.parse_input(sys.stdin.buffer)
        else:
            raise ValueError(f"Unknown format for Input: {type(Input)}")
        if self.verbose > 0:
            print(self.usecols, self.dim_cols)
            print(type(self.usecols), type(self.dim_cols))
        for data, dim in data_generator:
            self.write_chunk(data, dim)
        self.close()

    def catmz(self, Input=None, order=None):
        """
		Cat multiple .mz files into one .mz file.

		Parameters
		----------
		Input : str or list
			Either a str (including *, as input for glob, should be inside the
			double quotation marks if using fire) or a list.
		order : None, list or path
			If order=None, Input will be sorted using python sorted.
			If order is a list, tuple or array of basename.rstrip(.mz), sorted as order.
			If order is a file path (for example, chrom size path to order chroms
			or only use selected chroms) will be sorted as
			the 1st column of the input file path (without header, tab-separated).
			default is None

		Returns
		-------

		"""
        if type(Input) == str and '*' in Input:
            Input = glob.glob(Input)
        if type(Input) != list:
            raise ValueError("Input should be either a list of a string including *.")
        if order is None:
            Input = sorted(Input)
        else:
            D = {os.path.basename(inp)[:-3]: inp for inp in Input}
            if self.verbose > 0:
                print(D)
            if isinstance(order, str):
                order = pd.read_csv(os.path.abspath(os.path.expanduser(order)),
                                    sep='\t', header=None, usecols=[0])[0].tolist()
            if isinstance(order, (list, tuple, np.ndarray)):  # order is a list
                # Input=[str(i)+'.mz' for i in order]
                Input = [D[str(i)] for i in order]
            else:
                raise ValueError("input of order is not corrected !")
        if self.verbose > 0:
            print(Input)
        for file_path in Input:
            reader = Reader(file_path)
            data_size = reader.header['total_size'] - reader.header['header_size']
            self._handle.write(reader._handle.read(data_size))
            reader.close()
        self.write_ts_eof_close()

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

    def write_ts_eof_close(self):
        cur_offset = self._handle.tell()
        self._handle.seek(self._magic_size + 4)  # magic and version
        self._handle.write(struct.pack("<Q", cur_offset))  # real offset.
        self._handle.seek(cur_offset)
        self._handle.write(_bmz_eof)
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
        self.write_ts_eof_close()

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
if __name__ == "__main__":
    pass
