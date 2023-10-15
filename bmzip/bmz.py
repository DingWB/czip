#!/usr/bin/env python
import os,sys
import struct
import zlib
from builtins import open as _open
import numpy as np
import pandas as pd
import fire

_bmz_magic = b'BMZIP'
_block_magic = b"MB"
_chunk_magic = b"MC"
_BLOCK_MAX_LEN = 65535
_bmz_eof = b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BM\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
_version=1.0

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
	if block_start_offset < 0 or block_start_offset >= 281474976710656:
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
	if not magic or magic != _block_magic: # next chunk or EOF
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
		handle.seek(block_size-6, 1)
		data_len = struct.unpack("<H", handle.read(2))[0]
		return block_size, data_len
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
			handle = _open(Input, "rb")
		self._handle = handle
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
		self.header['magic']=magic
		self.header['version'] = struct.unpack("<f", f.read(4))[0]
		total_size=struct.unpack("<Q", f.read(8))[0]
		if total_size ==0:
			raise ValueError("File not completed !")
		self.header['total_size']=total_size
		format_len = struct.unpack("<B", f.read(1))[0]
		Formats=[]
		for i in range(format_len):
			n=struct.unpack("<B", f.read(1))[0]
			format = struct.unpack(f"<{n}s", f.read(n))[0].decode()
			Formats.append(format)
		self.header['Formats'] = Formats
		names=[]
		for i in range(format_len):
			n = struct.unpack("<B", f.read(1))[0]
			name = struct.unpack(f"<{n}s", f.read(n))[0].decode()
			names.append(name)
		self.header['names'] = names
		assert len(Formats) == len(names)
		tags=[]
		tag_len = struct.unpack("<B", f.read(1))[0]
		for i in range(tag_len):
			n = struct.unpack("<B", f.read(1))[0]
			tag = struct.unpack(f"<{n}s", f.read(n))[0].decode()
			tags.append(tag)
		self.header['tags'] = tags
		self.header['header_size'] = f.tell() #end of header, begin of 1st chunk
		self.fmt=''.join(Formats)
		self._unit_size=struct.calcsize(self.fmt)

	def print_header(self):
		print(self.header)

	def _load_chunk(self, start_offset=None, view=True):
		if start_offset is None: #this is another chunk, not the 1st one.
			start_offset = self._chunk_tail_offset
		if start_offset >= self.header['total_size']:
			return False
		self._handle.seek(start_offset)
		self._chunk_start_offset = start_offset # real offset on disk.
		magic = self._handle.read(2)
		if magic != _chunk_magic:
			raise StopIteration
		self._chunk_size = struct.unpack('<Q',self._handle.read(8))[0]
		self._chunk_tags=[]
		for t in self.header['tags']:
			n = struct.unpack("<B", self._handle.read(1))[0]
			tag = struct.unpack(f"<{n}s", self._handle.read(n))[0].decode()
			self._chunk_tags.append(tag)
		self._chunk_body_offset=self._handle.tell()
		# load chunk tail, jump CDATA part
		chunk_tail_start_offset = self._chunk_start_offset + self._chunk_size
		self._handle.seek(chunk_tail_start_offset)
		self._chunk_data_len=struct.unpack("<Q", self._handle.read(8))[0]
		self._chunk_nblocks=struct.unpack("<Q", self._handle.read(8))[0]
		self._chunk_block_offsets = []
		if view: # no need to load _chunk_block_offsets
			self._handle.seek(self._chunk_nblocks * 8,1)
			#_chunk_tail_offset = end position of this chunk
			# = start position of next chunk.
		else:
			# read block_offset
			for i in range(self._chunk_nblocks):
				block_offset=struct.unpack("<Q", self._handle.read(8))[0]
				self._chunk_block_offsets.append(block_offset)
		self._chunk_tail_offset=self._handle.tell()
		return True

	def get_chunks(self):
		r=self._load_chunk(self.header['header_size'],view=False)
		while r:
			yield [self._chunk_start_offset, self._chunk_size,
				   self._chunk_tags, self._chunk_data_len, self._chunk_tail_offset,
				   self._chunk_nblocks,self._chunk_block_offsets
				   ]
			r=self._load_chunk(view=False)

	def summary_blocks(self):
		r = self._load_chunk(self.header['header_size'], view=True)
		header=['chunk_start_offset','chunk_size','chunk_tags','chunk_data_len',
				'chunk_tail_offset','chunk_nblocks']+['block_start_offset',
			   'block_size', 'block_data_start', 'block_data_len']
		sys.stdout.write('\t'.join(header) + '\n')
		while r:
			self._handle.seek(self._chunk_body_offset)
			for block in SummaryBmzBlocks(self._handle):
				chunk_info=[self._chunk_start_offset, self._chunk_size,
				   self._chunk_tags, self._chunk_data_len, self._chunk_tail_offset,
				   self._chunk_nblocks]
				block=chunk_info+list(block)
				# print(block)
				#start_offset, block_length, data_start, data_len
				try:
					sys.stdout.write('\t'.join([str(v) for v in block])+'\n')
				except:
					sys.stdout.close()
			r = self._load_chunk(view=True)
		sys.stdout.close()

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
			block_size, self._buffer = _load_bmz_block(self._handle,True)
		except StopIteration: # EOF
			block_size = 0
			self._buffer = b""
		self._within_block_offset = 0
		self._block_raw_length = block_size

	def _print_cache(self):
		self._cached_data += self._buffer
		end_index = len(self._cached_data) - (len(self._cached_data) % self._unit_size)
		for result in struct.iter_unpack(f"<{self.fmt}", self._cached_data[:end_index]):
			line = '\t'.join([str(v) for v in result])
			try:
				sys.stdout.write(self.chunk_tag + line + '\n')
			except:
				sys.stdout.close()
		self._cached_data = self._cached_data[end_index:]

	def view(self, tag=None, header=True):
		"""
		View .bmz file.

		Parameters
		----------
		tag : str
			index of tags given to writer.write_chunk, separated by comma,
			default is None, tag will be shown in each row (such as sampleID
			and chrom)
		header : bool
			whether to print the header.

		Returns
		-------

		"""
		if type(tag)==int:
			tag=[tag]
		elif not tag is None:
			tag=[int(i) for i in tag.split(',')]

		if not tag is None:
			header_tag = "\t".join([self.header['tags'][t] for t in tag]) + '\t'
		else:
			header_tag = ''
		if header:
			line = "\t".join(self.header['names'])
			sys.stdout.write(header_tag + line + '\n')

		r = self._load_chunk(self.header['header_size'], view=True)
		# header_size is the start offset of 1st chunk
		while r: # chunk
			self._cached_data=b''
			# here, we have already read the chunk header and tail but
			# not the blocks, we need to jump to body start pos and read each block.
			self._load_block(start_offset = self._chunk_body_offset) #
			if not tag is None:
				self.chunk_tag="\t".join([self._chunk_tags[t] for t in tag])+'\t'
			else:
				self.chunk_tag=''
			while self._block_raw_length > 0:
				# deal with such case: unit_size is 10, but data(_buffer) size is 18,
				self._print_cache()
				self._load_block()
			r = self._load_chunk(view=True)
		sys.stdout.close()

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
			# (this avoids a function call since _load_block would do nothing)
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
		# assert virtual_offset == self.tell(), \
		#    "Did seek to %i (%i, %i), but tell says %i (%i, %i)" \
		#    % (virtual_offset, start_offset, within_block,
		#       self.tell(), self._block_start_offset,
		#       self._within_block_offset)
		return virtual_offset

	def read(self, size=-1):
		"""Read method for the BGZF module."""
		if size < 0:
			raise NotImplementedError("Don't be greedy, that could be massive!")

		result =  b""
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
				 names=['mc', 'cov'], tags=['chrom'],fileobj=None,
				 compresslevel=6,verbose=0):
		if Output and fileobj:
			raise ValueError("Supply either Output or fileobj, not both")
		if fileobj: # write to an existed openned file object
			if fileobj.read(0) != b"":
				raise ValueError("fileobj not opened in binary mode")
			handle = fileobj
		else:
			if not Output is None: #write output to a file
				if "w" not in mode.lower() and "a" not in mode.lower():
					raise ValueError(f"Must use write or append mode, not {mode!r}")
				handle = _open(Output, mode)
			else: # write to stdout buffer
				handle = sys.stdout.buffer
		self._handle = handle
		self._buffer = b""
		self._chunk_start_offset = None
		self._chunk_tags = None
		if type(Formats) == str and ',' not in Formats:
			self.Formats = [Formats]
		elif ',' in Formats:
			self.Formats = Formats.split(',')
		else:
			self.Formats = Formats
		self.fmts=''.join(list(Formats))
		self.compresslevel = compresslevel
		if type(names) == str and ',' not in names:
			self.names = [names]
		elif ',' in names:
			self.names = names.split(',')
		else:
			self.names = names
		if type(tags)==str and ',' not in tags:
			self.tags = [tags]
		elif ',' in tags:
			self.tags=tags.split(',')
		else:
			self.tags=tags
		self.magic_size = len(_bmz_magic)
		self.verbose = verbose
		self.write_header()

	def write_header(self):
		f=self._handle
		f.write(struct.pack("<5s", _bmz_magic))  # Identifier: 5 bytes
		f.write(struct.pack("<f", _version))  # 4 bytes, float
		f.write(struct.pack("<Q", 0)) # 8 bytes, total size, including magic.
		f.write(struct.pack("<B", len(self.Formats)))  # 1byte
		assert len(self.names) == len(self.Formats)
		for format in self.Formats:
			format_len = len(format)
			f.write(struct.pack("<B", format_len))  # length of each format, 1 byte
			f.write(struct.pack(f"<{format_len}s", bytes(format, 'utf-8')))
		for name in self.names:
			name_len = len(name)
			f.write(struct.pack("<B", name_len))  # length of each name, 1 byte
			f.write(struct.pack(f"<{name_len}s", bytes(name, 'utf-8')))
		f.write(struct.pack("<B", len(self.tags)))  # 1byte
		for tag in self.tags:
			tag_len = len(tag)
			f.write(struct.pack("<B", tag_len))  # length of each name, 1 byte
			f.write(struct.pack(f"<{tag_len}s", bytes(tag, 'utf-8')))

	def _write_block(self, block):
		if len(block) > _BLOCK_MAX_LEN: #65536 = 1 << 16 = 2**16
			raise ValueError(f"{len(block)} Block length > {_BLOCK_MAX_LEN}")
		c = zlib.compressobj(
			self.compresslevel, zlib.DEFLATED, -15, zlib.DEF_MEM_LEVEL, 0
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
		#block size: magic (2 btyes) + block_size (2bytes) + compressed data +
		# block_data_len (2 bytes)
		uncompressed_length = struct.pack("<H", len(block)) # 2 bytes
		data = _block_magic + bsize + compressed + uncompressed_length
		self._block_offsets.append(self._handle.tell())
		# _block_offsets are real position on disk, not a virtual offset
		self._handle.write(data)
		self._chunk_data_len += len(block)

	def _write_chunk_tail(self):
		# tail: chunk_data_len (Q,8bytes) + n_blocks (Q, 8bytes)
		# + list of block start offsets (n_blocks * 8 bytes)
		self._handle.write(struct.pack("<Q", self._chunk_data_len))
		self._handle.write(struct.pack("<Q", len(self._block_offsets)))
		for block_offset in self._block_offsets:
			self._handle.write(struct.pack("<Q", block_offset))

	def _chunk_finished(self):
		# current position is the end of chunk.
		cur_offset = self._handle.tell() # before chunk tail, not a virtual offset
		# go back and write the total chunk size
		self._handle.seek(self._chunk_start_offset + 2) #2 bytes for magic
		chunk_size = cur_offset - self._chunk_start_offset
		# chunk_size including the chunk_size itself, but not including chunk tail.
		self._handle.write(struct.pack("<Q", chunk_size))
		# go back the current position
		self._handle.seek(cur_offset) #not a virtual offset
		self._write_chunk_tail()

	def write_chunk(self, data, tags): # tags is a list.
		if isinstance(data, str):
			data = data.encode("latin-1")
		if self._chunk_tags != tags:
			# the first chunk or another new chunk
			if self._chunk_tags is not None:
				# this is another new chunk, current position is the end of chunk.
				self.flush() # finish_chunk in flush
			# else: #the first chunk, no chunk was writen previously.
			self._chunk_tags = tags
			self._chunk_start_offset = self._handle.tell()
			self._handle.write(_chunk_magic)
			# chunk total size place holder: 0
			self._handle.write(struct.pack("<Q", 0)) #8 bytes; including this chunk_size
			for tag in self._chunk_tags:
				tag_len = len(tag)
				self._handle.write(struct.pack("<B", tag_len))  # length of each name, 1 byte
				self._handle.write(struct.pack(f"<{tag_len}s", bytes(tag, 'utf-8')))
			self._chunk_data_len = 0
			self._block_offsets = []
			if self.verbose > 0:
				print("Writing chunk with tags: ",self._chunk_tags)

		data_len = len(data)
		if len(self._buffer) + data_len < _BLOCK_MAX_LEN:
			self._buffer += data
		else:
			self._buffer += data
			while len(self._buffer) >= _BLOCK_MAX_LEN:
				self._write_block(self._buffer[:_BLOCK_MAX_LEN])
				self._buffer = self._buffer[_BLOCK_MAX_LEN:]

	def parse_input(self,input_handle):
		if isinstance(input_handle,pd.DataFrame):
			# usecols and tag_cols should be in the columns of this dataframe.
			if self.chunksize is None:
				for tag, df1 in input_handle.groupby(self.tag_cols)[self.usecols]:
					if type(tag) != list:
						tag = [tag]
					yield df1.apply(lambda x: struct.pack(self.fmts, *x.tolist()),
									 axis=1).sum(),tag
			else:
				while input_handle.shape[0]>0:
					df=input_handle.iloc[:self.chunksize]
					for tag, df1 in df.groupby(self.tag_cols)[self.usecols]:
						if type(tag) != list:
							tag = [tag]
						yield df1.apply(lambda x: struct.pack(self.fmts, *x.tolist()),
										axis=1).sum(), tag
					input_handle=input_handle.iloc[self.chunksize:]
		else: # stdin or read from file.
			for df in pd.read_csv(input_handle,sep=self.sep,
								  usecols=self.tag_cols+self.usecols,
								  chunksize=self.chunksize,header=self.header,
								  skiprows=self.skiprows):
				for tag, df1 in df.groupby(self.tag_cols)[self.usecols]:
					if type(tag) != list:
						tag = [tag]
					yield df1.apply(lambda x: struct.pack(self.fmts, *x.tolist()),
									 axis=1).sum(),tag

	def pack(self, Input, usecols=[1,4,5],tag_cols=[0],
			 sep='\t',chunksize=5000,header=None, skiprows=0):
		self.sep = sep
		if type(usecols)==int or (type(usecols)==str and ',' not in usecols):
			self.usecols = [int(usecols)]
		elif ',' in usecols:
			self.usecols = [int(i) for i in usecols.split(',')]
		else:
			self.usecols = [int(i) for i in usecols]
		if type(tag_cols)==int or (type(tag_cols)==str and ',' not in tag_cols):
			self.tag_cols = [int(tag_cols)]
		elif ',' in tag_cols:
			self.tag_cols = [int(i) for i in tag_cols.split(',')]
		else:
			self.tag_cols = [int(i) for i in tag_cols]
		self.chunksize = chunksize
		self.header = header
		self.skiprows = skiprows
		if isinstance(Input,(list,tuple,np.ndarray)):
			Input=pd.DataFrame(Input)
		if Input=='stdin' or Input=='-':
			data_generator=self.parse_input(sys.stdin.buffer)
		elif isinstance(Input, str):
			input_path = os.path.abspath(os.path.expanduser(Input))
			data_generator = self.parse_input(input_path)
		else:  # Input is a dataframe
			data_generator = self.parse_input(Input)

		for data,tag in data_generator:
			self.write_chunk(data,tag)
		self.close()

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

	def save(self):
		if self._buffer:
			self.flush()
		# write total size in the file header, just after magic (5 bytes)
		cur_offset = self._handle.tell()
		self._handle.seek(9) #5 is magic and version (5+4=9 bytes)
		self._handle.write(struct.pack("<Q", cur_offset)) # real offset.
		self._handle.seek(cur_offset)

	def close(self):
		"""Flush data, write 28 bytes BGZF EOF marker, and close BGZF file.

		samtools will look for a magic EOF marker, just a 28 byte empty BGZF
		block, and if it is missing warns the BAM file may be truncated. In
		addition to samtools writing this block, so too does bgzip - so this
		implementation does too.
		"""
		self.save()
		self._handle.write(_bmz_eof)
		self._handle.flush()
		self._handle.close()

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