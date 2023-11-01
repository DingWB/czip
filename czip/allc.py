#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:34:38 2020

@author: DingWB
"""
import itertools
import sys
import os, sys
import struct
import pandas as pd
import pysam
import multiprocessing
from Bio import SeqIO
import numpy as np
from collections import defaultdict
import math
from .cz import (Reader, Writer, get_dtfuncs,
                 _BLOCK_MAX_LEN, _chunk_magic)


# ==========================================================
def WriteC(record, outdir, chunksize=5000):
    # cdef int i, N
    # cdef char* chrom, base, context, strand
    chrom = record.id
    outfile = os.path.join(outdir, chrom + ".cz")
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return None
    print(chrom)
    writer = Writer(outfile, Formats=['Q', 'c', '3s'],
                    Columns=['pos', 'strand', 'context'],
                    Dimensions=['chrom'])
    dtfuncs = get_dtfuncs(writer.Formats)
    N = record.seq.__len__()
    data = b''
    for i in range(N):  # 0-based
        base = record.seq[i:i + 1].upper()
        if base.__str__() == 'C':  # forward strand
            context = record.seq[i: i + 3].upper().__str__()  # pos, left l1 base pair and right l2 base pair
            strand = '+'
        elif base.reverse_complement().__str__() == 'C':  # reverse strand
            context = record.seq[i - 2:i + 1].reverse_complement().upper().__str__()
            strand = '-'
        else:
            continue
        # f.write(f"{chrom}\t{i}\t{i + 1}\t{context}\t{strand}\n")
        values = [func(v) for v, func in zip([i + 1, strand, context], dtfuncs)]
        data += struct.pack(writer.fmts, *values)
        # position is 0-based (start) 1-based (end position, i+1)
        if i % chunksize == 0 and len(data) > 0:
            writer.write_chunk(data, [chrom])
            data = b''
    if len(data) > 0:
        writer.write_chunk(data, [chrom])
    writer.close()


# ==========================================================
class AllC:
    def __init__(self, Genome=None, Output="hg38_allc.cz",
                 pattern="C", n_jobs=12, keep_temp=False):
        """
        Extract position of specific pattern in the reference genome, for example C.
            Example: python ~/Scripts/python/tbmate.py AllC -g ~/genome/hg38/hg38.fa --n_jobs 10 run
            Or call within python: ac=AllC(genome="/gale/netapp/home2/wding/genome/hg38/hg38.fa")
        Parameters
        ----------
        genome: path
            reference genome.
        out: path
            path for output
        pattern: str
            pattern [C]
        n_jobs: int
            number of CPU used for Pool.
        """
        self.genome=os.path.abspath(os.path.expanduser(Genome))
        self.Output=os.path.abspath(os.path.expanduser(Output))
        self.outdir=self.Output+'.tmp'
        if not os.path.exists(self.outdir):
            os.mkdir(self.outdir)
        self.pattern=pattern
        self.records = SeqIO.parse(self.genome, "fasta")
        self.n_jobs = n_jobs if not n_jobs is None else os.cpu_count()
        self.keep_temp = keep_temp
        if pattern=='C':
            self.func=WriteC

    def writePattern(self):
        pool = multiprocessing.Pool(self.n_jobs)
        jobs = []
        for record in self.records:
            job = pool.apply_async(self.func, (record,self.outdir))
            jobs.append(job)
        for job in jobs:
            job.get()
        pool.close()
        pool.join()

    def merge(self):
        writer = Writer(Output=self.Output, Formats=['Q', 'c', '3s'],
                        Columns=['pos', 'strand', 'context'],
                        Dimensions=['chrom'], message=self.genome)
        writer.catcz(Input=f"{self.outdir}/*.cz")

    def run(self):
        self.writePattern()
        self.merge()
        if not self.keep_temp:
            os.system(f"rm -rf {self.outdir}")


def bed2cz(input, outfile, reference=None, missing_value=[0, 0],
           Formats=['H', 'H'], Columns=['mc', 'cov'], Dimensions=['chrom'],
           usecols=[4, 5], pr=0, pa=1, sep='\t', Path_to_chrom=None,
           chunksize=2000):
    """
    convert allc.tsv.gz to .cz file.

    Parameters
    ----------
    input : path
        path to allc.tsv.gz, should has .tbi index.
    outfile : path
        output .cz file
    reference : path
        path to reference coordinates.
    Formats: list
        When reference is provided, we only need to pack mc and cov,
        ['H', 'H'] is suggested (H is unsigned short integer, only 2 bytes),
        if reference is not provided, we also need to pack position (Q is
        recommanded), in this case, Formats should be ['Q','H','H'].
    Columns: list
        Columns names, in default is ['mc','cov'] (reference is provided), if no
        referene provided, one should use ['pos','mc','cov'].
    Dimensions: list
        Dimensions passed to czip.Writer, dimension name, for allc file, dimension
        is chrom.
    usecols: list
        default is [4, 5], for a typical .allc.tsv.gz, if no reference is provided,
        the columns to be packed should be [1,4,5] (pos, mv and cov).
        If reference is provided, then we only need to pack [4,5] (mc and cov).
    pr: int
        index of position column in reference .cz header columns [0]
    pa: int
        index of position column in input input or bed column.
    chunksize : int
        default is 5000
    Path_to_chrom : path
        path to chrom_size path or similar file containing chromosomes order,
        the first columns should be chromosomes, tab separated and no header.

    Returns
    -------

    """
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    allc_path = os.path.abspath(os.path.expanduser(input))
    if not os.path.exists(allc_path + '.tbi'):
        raise ValueError(f"index file .tbi not existed, please create index first.")
    print(allc_path)
    tbi = pysam.TabixFile(allc_path)
    contigs = tbi.contigs
    if not Path_to_chrom is None:
        Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
        df = pd.read_csv(Path_to_chrom, sep='\t', header=None, usecols=[0])
        chroms = df.iloc[:, 0].tolist()
        all_chroms = [c for c in chroms if c in contigs]
    else:
        all_chroms = contigs
    if not reference is None:
        reference = os.path.abspath(os.path.expanduser(reference))
        message = os.path.basename(reference)
    else:
        message = ''
        # formats, columns, dimensions = ['Q', 'H', 'H'], ['pos', 'mc', 'cov'], ['chrom']
        # usecols = [1, 4, 5]
    writer = Writer(outfile, Formats=Formats, Columns=Columns,
                    Dimensions=Dimensions, message=message)
    dtfuncs = get_dtfuncs(Formats, tobytes=False)

    if not reference is None:
        ref_reader = Reader(reference)
        na_value_bytes = struct.pack(f"<{writer.fmts}", *missing_value)
        for chrom in all_chroms:
            # print(chrom,'\t'*4,end='\r')
            # method 1
            ref_positions = ref_reader.__fetch__(tuple([chrom]), s=pr, e=pr + 1)
            records = tbi.fetch(chrom)
            data, i = b'', 0
            row_query = next(records).rstrip('\n').split(sep)
            row_query_pos = int(row_query[pa])
            for ref_pos in ref_positions:
                if ref_pos[0] < row_query_pos:
                    data += na_value_bytes
                    i += 1
                else:
                    if ref_pos[0] == row_query_pos:  # match
                        data += struct.pack(f"<{writer.fmts}",
                                            *[func(row_query[i]) for i, func in
                                              zip(usecols, dtfuncs)
                                              ])
                        i += 1
                    try:
                        row_query = next(records).rstrip('\n').split(sep)
                        row_query_pos = int(row_query[pa])
                    except:
                        break
                if i > chunksize:
                    writer.write_chunk(data, [chrom])
                    data, i = b'', 0
            # when break: write all the rest reference positions (na_value_bytes)
            for ref_pos in ref_positions:
                data += na_value_bytes
                i += 1
                if i > chunksize:
                    writer.write_chunk(data, [chrom])
                    data, i = b'', 0
            if len(data) > 0:
                writer.write_chunk(data, [chrom])

            # methd 2, similar speed, but larger memory requirement.
            # ref_positions = np.array([pos[0] for pos in ref_reader.__fetch__(tuple([chrom]), s=pr, e=pr + 1)])
            # records = [line.rstrip('\n').split(sep) for line in tbi.fetch(chrom)]
            # query_positions = np.array([int(record[pa]) for record in records])
            # indices = np.where(np.in1d(ref_positions, query_positions))[0]
            # # indices is the indice where element of query_positions in ref_positions
            # indice_start = 0
            # # sum=0
            # for indice, record in zip(indices, records):
            #     for i in range((indice - indice_start) // chunksize):
            #         writer.write_chunk(na_value_bytes * chunksize, [chrom])
            #         # sum+=chunksize
            #     i = (indice - indice_start) % chunksize
            #     data = b''
            #     if i > 0:
            #         data += na_value_bytes * i
            #         # sum += i
            #     data += struct.pack(f"<{writer.fmts}", *[func(record[i]) for i, func in
            #                                              zip(usecols, dtfuncs)])
            #     writer.write_chunk(data, [chrom])
            #     # sum += 1
            #     indice_start = indice + 1
            # indice = len(ref_positions)
            # for i in range((indice - indice_start) // chunksize):
            #     writer.write_chunk(na_value_bytes * chunksize, [chrom])
            #     # sum += chunksize
            # i = (indice - indice_start) % chunksize
            # if i > 0:
            #     writer.write_chunk(na_value_bytes * i, [chrom])
            #     # sum += i
        ref_reader.close()
    else:
        for chrom in all_chroms:
            # print(chrom)
            # records = tbi.fetch(chrom)
            # data, i = b'', 0
            # while True:
            #     try:
            #         values = records.__next__().rstrip('\n').split(sep)
            #     except:
            #         break
            #     if i >= chunksize:  # dims are the same, but reach chunksize
            #         writer.write_chunk(data, [chrom])
            #         data, i = b'', 0
            #     values = [func(values[i]) for i, func in zip(usecols, dtfuncs)]
            #     data += struct.pack(f"<{writer.fmts}", *values)
            #     i += 1
            # if len(data) > 0:
            #     writer.write_chunk(data, [chrom])

            # method 2
            for line in tbi.fetch(chrom):
                values = line.rstrip('\n').split(sep)
                data = struct.pack(f"<{writer.fmts}",
                                   *[func(values[i]) for i, func in zip(usecols, dtfuncs)])
                writer.write_chunk(data, [chrom])
    writer.close()
    tbi.close()


# ==========================================================
def _isCG(record):
    return record[2][:2] == b'CG'


# ==========================================================
def _isForwardCG(record):
    return record[2][:2] == b'CG' and record[1] == b'+'


# ==========================================================
def _isCH(record):
    return not _isCG(record)


# ==========================================================
def generate_ssi(Input, output=None, pattern="CGN"):
    """
    Generate ssi (subset index) for a given input .cz

    Parameters
    ----------
    Input : .cz
    output : .ssi
    pattern : CGN, CHN, +CGN, -CGN

    Returns
    -------

    """
    if pattern == 'CGN':
        judge_func = _isCG
    elif pattern == 'CHN':  # CH
        judge_func = _isCH
    elif pattern == '+CGN':
        judge_func = _isForwardCG
    else:
        raise ValueError("Currently, only CGN, CHN, +CGN supported")
    if output is None:
        output = Input + '.' + pattern + '.ssi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(Input)
    reader.category_ssi(output=output, formats=['I'], columns=['ID'],
                        dimensions=['chrom'], match_func=judge_func,
                        chunksize=2000)
    reader.close()


def merge_cz_worker(outfile_cat, outdir, chrom, dims, formats,
                    block_idx_start, batch_nblock, chunksize=5000):
    # print(chrom, "started",end='\r')
    outname = os.path.join(outdir, chrom + f'.{block_idx_start}.cz')
    reader1 = Reader(outfile_cat)
    data = None
    for dim in dims:  # each dim is a file of the same chrom
        r = reader1._load_chunk(reader1.dim2chunk_start[dim], jump=False)
        block_start_offset = reader1._chunk_block_1st_record_virtual_offsets[
                                 block_idx_start] >> 16
        buffer = b''
        for i in range(batch_nblock):
            reader1._load_block(start_offset=block_start_offset)
            buffer += reader1._buffer
            block_start_offset = None
        values = np.array([result[:2] for result in struct.iter_unpack(
            f"<{reader1.fmts}", buffer)])  # values for batch_nblock (23) blocks
        if data is None:
            data = values.copy()
        else:
            data += values
        # q.put(tuple([dim, chunk_id, data]))
    writer1 = Writer(outname, Formats=formats,
                     Columns=reader1.header['Columns'],
                     Dimensions=reader1.header['Dimensions'][:1],
                     message=outfile_cat)
    byte_data, i = b'', 0
    for values in data.tolist():
        byte_data += struct.pack(f"<{writer1.fmts}",
                                 *values)
        i += 1
        if i > chunksize:
            writer1.write_chunk(byte_data, [chrom])
            byte_data, i = b'', 0
    if len(byte_data) > 0:
        writer1.write_chunk(byte_data, [chrom])
    writer1.close()
    reader1.close()
    print(chrom, block_idx_start, "done")
    return


def merge_cz(indir=None, cz_paths=None, outfile="merged.cz", n_jobs=12, formats=['I', 'I'],
             Path_to_chrom="~/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt",
             keep_cat=False, batchsize=10):
    outfile = os.path.abspath(os.path.expanduser(outfile))
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    if cz_paths is None:
        cz_paths = os.listdir(indir)
    reader = Reader(os.path.join(indir, cz_paths[0]))
    header = reader.header
    reader.close()
    outfile_cat = outfile + '.cat.cz'
    writer = Writer(Output=outfile_cat, Formats=header['Formats'],
                    Columns=header['Columns'], Dimensions=header['Dimensions'],
                    message="catcz")
    writer.catcz(Input=[os.path.join(indir, cz_path) for cz_path in cz_paths],
                 add_dim=True)
    Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
    df = pd.read_csv(Path_to_chrom, sep='\t', header=None, usecols=[0])
    chroms = df.iloc[:, 0].tolist()

    reader = Reader(outfile_cat)
    chrom_col = reader.header['Dimensions'][0]
    chunk_info = reader.chunk_info
    reader.close()
    chrom_nblocks = chunk_info.reset_index().loc[:, [chrom_col, 'chunk_nblocks']
                    ].drop_duplicates().set_index(chrom_col).chunk_nblocks.to_dict()
    # chunk_info.set_index(chrom_col)
    unit_nblock = int(writer._unit_size / (math.gcd(writer._unit_size, _BLOCK_MAX_LEN)))
    nunit_perbatch = int(np.ceil((chunk_info.chunk_nblocks.max() / batchsize
                                  ) / unit_nblock))
    batch_nblock = nunit_perbatch * unit_nblock  # how many block for each batch
    # manager = multiprocessing.Manager()
    # queue1 = manager.Queue()
    pool = multiprocessing.Pool(n_jobs)
    # watcher = pool.apply_async(czs_merger, (outfile, formats, columns,
    #                                         dimensions, message, n_batch, queue1))
    jobs = []
    outdir = outfile + '.tmp'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for chrom in chroms:
        dims = chunk_info.loc[chunk_info[chrom_col] == chrom].index.tolist()
        if len(dims) == 0:
            continue
        block_idx_start = 0
        while block_idx_start < chrom_nblocks[chrom]:
            job = pool.apply_async(merge_cz_worker,
                                   (outfile_cat, outdir, chrom, dims, formats,
                                    block_idx_start, batch_nblock))
            jobs.append(job)
            block_idx_start += batch_nblock
    for job in jobs:
        r = job.get()
    # queue1.put('kill')
    pool.close()
    pool.join()
    # manager._process.terminate()
    # manager.shutdown()
    if not keep_cat:
        os.remove(outfile_cat)
    # First, merge different batch for each chrom
    for chrom in chroms:
        # print(chrom,end='\t')
        outname = os.path.join(outdir, chrom + '.cz')
        writer = Writer(Output=outname, Formats=formats,
                        Columns=header['Columns'], Dimensions=header['Dimensions'],
                        message=outfile_cat)
        writer._chunk_start_offset = writer._handle.tell()
        writer._handle.write(_chunk_magic)
        # chunk total size place holder: 0
        writer._handle.write(struct.pack("<Q", 0))  # 8 bytes; including this chunk_size
        writer._chunk_data_len = 0
        writer._block_1st_record_virtual_offsets = []
        writer._chunk_dims = [chrom]
        block_idx_start = 0
        infile = os.path.join(outdir, chrom + f'.{block_idx_start}.cz')
        while os.path.exists(infile):
            reader = Reader(infile)
            reader._load_chunk(reader.header['header_size'])
            block_start_offset = reader._chunk_start_offset + 10
            writer._buffer = b''
            for i in range(reader._chunk_nblocks):
                reader._load_block(start_offset=block_start_offset)  #
                if len(writer._buffer) + len(reader._buffer) < _BLOCK_MAX_LEN:
                    writer._buffer += reader._buffer
                else:
                    writer._buffer += reader._buffer
                    while len(writer._buffer) >= _BLOCK_MAX_LEN:
                        writer._write_block(writer._buffer[:_BLOCK_MAX_LEN])
                        writer._buffer = writer._buffer[_BLOCK_MAX_LEN:]
                block_start_offset = None
            block_idx_start += batch_nblock
            infile = os.path.join(outdir, chrom + f'.{block_idx_start}.cz')
            reader.close()
        # write chunk tail
        writer.close()

    # Second, merge chromosomes to outfile
    writer = Writer(Output=outfile, Formats=formats,
                    Columns=header['Columns'], Dimensions=header['Dimensions'],
                    message="merged")
    writer.catcz(Input=[f"{outdir}/{chrom}.cz" for chrom in chroms])
    os.system(f"rm -rf {outdir}")


def merge_cell_type(indir=None, cell_table=None, outdir=None,
                    n_jobs=64, Path_to_chrom=None, ext='.CGN.merged.cz'):
    indir = os.path.abspath(os.path.expanduser(indir))
    outdir = os.path.abspath(os.path.expanduser(outdir))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
    df_ct = pd.read_csv(cell_table, sep='\t', header=None, names=['cell', 'ct'])
    for ct in df_ct.ct.unique():
        outfile = os.path.join(outdir, ct + '.cz')
        if os.path.exists(outfile):
            print(f"{outfile} existed.")
            continue
        print(ct)
        snames = df_ct.loc[df_ct.ct == ct, 'cell'].tolist()
        cz_paths = [os.path.join(indir, sname + ext) for sname in snames]
        merge_cz(indir=indir, cz_paths=cz_paths,
                 outfile=outfile, n_jobs=n_jobs, Path_to_chrom=Path_to_chrom)
# ==========================================================
def extractCG(input=None, outfile=None, ssi=None, chunksize=5000,
              merge_strand=True):
    """
    Extract CG context from .cz file

    Parameters
    ----------
    cz_path :path
    outfile :path
    ssi : path
        ssi should be ssi to mm10_with_chrL.allc.cz.CGN.ssi, not forward
        strand ssi, but after merge (if merge_strand is True), forward ssi
        mm10_with_chrL.allc.cz.+CGN.ssi should be used to generate
         reference, one can
        run: czip extract -m mm10_with_chrL.allc.cz
        -o mm10_with_chrL.allCG.forward.cz
        -b mm10_with_chrL.allc.cz.+CGN.ssi and use
        mm10_with_chrL.allCG.forward.cz as new reference.
    chunksize :int
    merge_strand: bool
        after merging, only forward strand would be kept, reverse strand values
        would be added to the corresponding forward strand.

    Returns
    -------

    """
    cz_path = os.path.abspath(os.path.expanduser(input))
    ssi_path = os.path.abspath(os.path.expanduser(ssi))
    ssi_reader = Reader(ssi_path)
    reader = Reader(cz_path)
    writer = Writer(outfile, Formats=reader.header['Formats'],
                    Columns=reader.header['Columns'],
                    Dimensions=reader.header['Dimensions'],
                    message=ssi_path)
    dtfuncs = get_dtfuncs(writer.Formats)
    for dim in reader.dim2chunk_start.keys():
        # print(dim)
        IDs = ssi_reader.get_ids_from_ssi(dim)
        if len(IDs.shape) != 1:
            raise ValueError("Only support 1D ssi now!")
        records = reader._getRecordsByIds(dim, IDs)
        data, count = b'', 0
        # for CG, if pos is forward (+), then pos+1 is reverse strand (-)
        if merge_strand:
            for i, record in enumerate(records):  # unpacked bytes
                if i % 2 == 0:
                    v0 = struct.unpack(f"<{reader.fmts}", record)
                else:
                    v1 = struct.unpack(f"<{reader.fmts}", record)
                    values = [r1 + r2 for r1, r2 in zip(v0, v1)]
                    data += struct.pack(writer.fmts,
                                        *[func(v) for v, func in zip(values, dtfuncs)])
                    count += 1
                if count > chunksize:
                    writer.write_chunk(data, dim)
                    data, count = b'', 0
        else:
            for record in records:  # unpacked bytes
                data += record
                count += 1
                if count > chunksize:
                    writer.write_chunk(data, dim)
                    data, count = b'', 0
        if len(data) > 0:
            writer.write_chunk(data, dim)
    writer.close()
    reader.close()
    ssi_reader.close()
# ==========================================================
if __name__ == "__main__":
    import fire
    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire()
