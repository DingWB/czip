#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:34:38 2020

@author: DingWB
"""
import itertools
import sys
import os
import struct
import pandas as pd
import pysam
import multiprocessing
from Bio import SeqIO
from collections import defaultdict
import random
# import pyximport
# pyximport.install(pyimport=True) #pyximport.install(pyimport=True)
from .utils import WriteC
from .bmz import Reader, Writer, get_dtfuncs
import numba
# ==========================================================
class AllC:
    def __init__(self, Genome=None, Output="hg38_allc.mz",
                 pattern="C", n_jobs=12):
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
        writer.catmz(Input=f"{self.outdir}/*.mz")

    def run(self):
        self.writePattern()
        self.merge()
        os.system(f"rm -rf {self.outdir}")


# ==========================================================
def allc2mz_worker_with_ref(allc_path, outdir, reference, chrom, all_chroms,
                            formats=['H', 'H'], columns=['mc', 'cov'],
                            dimensions=['chrom'], usecols=[4, 5],
                            missing_value=[0, 0], chunksize=2000,
                            pr=0, pa=1, sep='\t', q=None):
    """
    Pack allc to .mz file for one chrom, allc_path must has already been indexed
    using tabix.

    Parameters
    ----------
    allc_path :path
        path to allc file.
    outfile : path
        temporary .mz file, should be merged together after the Pool finished.
    reference: path
        path to reference .mz file. By providing reference, the pos in allc file will
        not be saved into .mz file, instead, we will align the position to the allc
        coordinates from the reference allc coordinates.
    chrom : chrom
        chrom, will be used as dimension to write chunk.
    verbose : int
        whether print debug information
    formats : list
        For allc file, default is ['Q','H', 'H'] for pos, mc, cov respectively.
    columns : list
        columsn names, default is [pos, mc, cov]
    dimensions : list
        For allc file, default dimension is ['chrom'].
    missing_value: list
        when the methylation record for a C is missing in allc file, then missing_value
        would be written into .mz file, default is [0,0] for mc and cov respectively.
    chunksize : int
        default is 5000
    usecols : list
        list of column index in the input file columns, in default, for allc file,
        usecols = [1,4,5], means [pos, mc, cov], 0-based.
    sep : str
        separator for the input file ['\t']

    Returns
    -------

    """
    tbi = pysam.TabixFile(allc_path)
    records = tbi.fetch(chrom)
    dtfuncs = get_dtfuncs(formats, tobytes=False)
    outfile = os.path.join(outdir, chrom + '.mz')
    writer = Writer(outfile, Formats=formats,
                    Columns=columns, Dimensions=dimensions)
    ref_reader = Reader(reference)
    ref_positions = ref_reader.__fetch__(tuple([chrom]), s=pr, e=pr + 1)
    data = b''
    na_value_bytes = struct.pack(f"<{writer.fmts}", *missing_value)
    i = 0
    while True:
        try:
            row_query = records.__next__().rstrip('\n').split(sep)
        except:
            break
        ref_pos = ref_positions.__next__()[0]
        while ref_pos < int(row_query[pa]):
            data += na_value_bytes
            i += 1
            try:
                ref_pos = ref_positions.__next__()[0]
            except:
                break
        if ref_pos == int(row_query[pa]):  # match
            values = [func(row_query[i]) for i, func in zip(usecols, dtfuncs)]
            data += struct.pack(f"<{writer.fmts}", *values)
            i += 1

        if i >= chunksize:
            writer.write_chunk(data, [chrom])
            data = b''
            i = 0
    if len(data) > 0:
        writer.write_chunk(data, [chrom])
    writer.close()
    ref_reader.close()
    tbi.close()
    # print(f"{outfile} done.", "\t" * 4, end='\r')
    if not q is None:
        q.put(tuple([outdir, chrom, all_chroms]))
    return outdir, chrom, all_chroms


# ==========================================================
def allc2mz_worker_without_ref(allc_path, outdir, chrom, all_chroms,
                               formats=['Q', 'H', 'H'], columns=['pos', 'mc', 'cov'],
                               dimensions=['chrom'], usecols=[1, 4, 5],
                               chunksize=2000, sep='\t', q=None):
    """
    Pack allc to .mz file for one chrom, allc_path must has already been indexed
    using tabix.

    Parameters
    ----------
    allc_path :path
        path to allc file.
    outfile : path
        temporary .mz file, should be merged together after the Pool finished.
    chrom : chrom
        chrom, will be used as dimension to write chunk.
    verbose : int
        whether print debug information
    formats : list
        For allc file, default is ['Q','H', 'H'] for pos, mc, cov respectively.
    columns : list
        columsn names, default is [pos, mc, cov]
    dimensions : list
        For allc file, default dimension is ['chrom'].
    chunksize : int
        default is 5000
    usecols : list
        list of column index in the input file columns, in default, for allc file,
        usecols = [1,4,5], means [pos, mc, cov], 0-based.
    sep : str
        separator for the input file ['\t']

    Returns
    -------

    """
    dtfuncs = get_dtfuncs(formats, tobytes=False)
    outfile = os.path.join(outdir, chrom + '.mz')
    tbi = pysam.TabixFile(allc_path)
    records = tbi.fetch(chrom)
    writer = Writer(outfile, Formats=formats,
                    Columns=columns, Dimensions=dimensions)
    data, i = b'', 0
    while True:
        try:
            values = records.__next__().rstrip('\n').split(sep)
        except:
            break
        if i >= chunksize:  # dims are the same, but reach chunksize
            writer.write_chunk(data, [chrom])
            data, i = b'', 0
        values = [func(values[i]) for i, func in zip(usecols, dtfuncs)]
        data += struct.pack(f"<{writer.fmts}", *values)
        i += 1
    if len(data) > 0:
        writer.write_chunk(data, [chrom])
    writer.close()
    tbi.close()
    # print(f"{outfile} done.", "\t" * 4, end='\r')
    if not q is None:
        q.put(tuple([outdir, chrom, all_chroms]))
    return outdir, chrom, all_chroms


def mzs_merger(formats, columns, dimensions, message, q):
    finished_jobs = {}
    while 1:
        result = q.get()  # if q is empty, it will wait.
        if result == 'kill':
            print("Done!")
            break
        outdir, chrom, all_chroms = result
        if outdir not in finished_jobs:
            finished_jobs[outdir] = [chrom]
        else:
            finished_jobs[outdir].append(chrom)
        if len(finished_jobs[outdir]) < len(all_chroms):
            continue  # wait
        outfile = outdir[:-4]  # .tmp
        writer = Writer(Output=outfile, Formats=formats,
                        Columns=columns, Dimensions=dimensions,
                        message=message)
        writer.catmz(Input=f"{outdir}/*", dim_order=all_chroms)
        os.system(f"rm -rf {outdir}")
        print(os.path.basename(outfile))
        del finished_jobs[outdir]  # finished one allc, remove from monitoring
    return

# ==========================================================
@numba.njit
def allc2mz(allc_path, outfile, reference=None, missing_value=[0, 0],
            chunksize=2000, pr=0, pa=1, sep='\t',
            Path_to_chrom=None):
    """
    convert allc.tsv.gz to .mz file.

    Parameters
    ----------
    allc_path : path
        path to allc.tsv.gz, should has .tbi index.
    outfile : path
        output .mz file
    reference : path
        path to reference coordinates.
    chunksize : int
        default is 5000
    Path_to_chrom : path
        path to chrom_size path or similar file containing chromosomes order,
        the first columns should be chromosomes, tab separated and no header.

    Returns
    -------

    """
    global row_query, records, fmts, dtfuncs, usecols
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    allc_path = os.path.abspath(os.path.expanduser(allc_path))
    if not os.path.exists(allc_path + '.tbi'):
        raise ValueError(f"index file .tbi not existed, please create index first.")
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
        formats, columns, dimensions = ['H', 'H'], ['mc', 'cov'], ['chrom']
        usecols = [4, 5]
    else:
        message = ''
        formats, columns, dimensions = ['Q', 'H', 'H'], ['pos', 'mc', 'cov'], ['chrom']
        usecols = [1, 4, 5]
    writer = Writer(outfile, Formats=formats, Columns=columns,
                    Dimensions=dimensions, message=message)
    dtfuncs = get_dtfuncs(formats, tobytes=False)

    def row2byte():
        global row_query, records, fmts, dtfuncs, usecols
        values = [func(row_query[i]) for i, func in zip(usecols, dtfuncs)]
        try:
            row_query = records.__next__().rstrip('\n').split('\t')
            return struct.pack(f"<{fmts}", *values)
        except:
            return False

    if not reference is None:
        ref_reader = Reader(reference)
        na_value_bytes = struct.pack(f"<{writer.fmts}", *missing_value)
        for chrom in all_chroms:
            print(chrom + '\t' * 10, end='\r')
            records = tbi.fetch(chrom)

            # method 2: better
            # row_query = records.__next__().rstrip('\n').split(sep)
            # data=b''.join(
            #     list(
            #         itertools.takewhile(lambda x:x!=False, [
            #             na_value_bytes if ref_pos[0] < int(row_query[pa])
            #             else row2byte()
            #             for ref_pos in ref_reader.__fetch__(
            #                 tuple([chrom]), s=pr, e=pr + 1
            #             )
            #         ])
            #     )
            # )
            # writer.write_chunk(data, [chrom])

            # method 1
            # ref_positions = ref_reader.__fetch__(tuple([chrom]), s=pr, e=pr + 1)
            # data = b''
            # i = 0
            # while True:
            #     try:
            #         row_query = records.__next__().rstrip('\n').split(sep)
            #     except:
            #         break
            #     ref_pos = ref_positions.__next__()[0]
            #     while ref_pos < int(row_query[pa]):
            #         data += na_value_bytes
            #         i += 1
            #         try:
            #             ref_pos = ref_positions.__next__()[0]
            #         except:
            #             break
            #     if ref_pos == int(row_query[pa]):  # match
            #         values = [func(row_query[i]) for i, func in zip(usecols, dtfuncs)]
            #         data += struct.pack(f"<{writer.fmts}", *values)
            #         i += 1
            #
            #     if i >= chunksize:
            #         writer.write_chunk(data, [chrom])
            #         data = b''
            #         i = 0
            # if len(data) > 0:
            #     writer.write_chunk(data, [chrom])

            # method 3
            ref_positions = ref_reader.__fetch__(tuple([chrom]), s=pr, e=pr + 1)
            data = b''
            while True:
                try:
                    row_query = records.__next__().rstrip('\n').split(sep)
                except:
                    break
                ref_pos = ref_positions.__next__()[0]
                while ref_pos < int(row_query[pa]):
                    data += na_value_bytes
                    try:
                        ref_pos = ref_positions.__next__()[0]
                    except:
                        break
                if ref_pos == int(row_query[pa]):  # match
                    values = [func(row_query[i]) for i, func in zip(usecols, dtfuncs)]
                    data += struct.pack(f"<{writer.fmts}", *values)
            writer.write_chunk(data, [chrom])
        ref_reader.close()
    else:
        for chrom in all_chroms:
            records = tbi.fetch(chrom)
            data, i = b'', 0
            while True:
                try:
                    values = records.__next__().rstrip('\n').split(sep)
                except:
                    break
                if i >= chunksize:  # dims are the same, but reach chunksize
                    writer.write_chunk(data, [chrom])
                    data, i = b'', 0
                values = [func(values[i]) for i, func in zip(usecols, dtfuncs)]
                data += struct.pack(f"<{writer.fmts}", *values)
                i += 1
            if len(data) > 0:
                writer.write_chunk(data, [chrom])
    writer.close()
    tbi.close()


# ==========================================================
def allc2mz_mp(allc_paths, outdir, reference=None, missing_value=[0, 0],
               chunksize=2000, pr=0, pa=1, n_jobs=8, sep='\t',
               Path_to_chrom=None, ext='.allc.tsv.gz'):
    """
    convert allc.tsv.gz to .mz file.

    Parameters
    ----------
    allc_paths : path
        path to allc.tsv.gz, should has .tbi index.
    outdir : path
        outdir to .mz file
    reference : path
        path to reference coordinates.
    n_jobs : int
        number of cpu to process each allc file
    chunksize : int
        default is 5000
    verbose : int
    Path_to_chrom : path
        path to chrom_size path or similar file containing chromosomes order,
        the first columns should be chromosomes, tab separated and no header.

    Returns
    -------

    """
    if isinstance(allc_paths, str):
        input_path = os.path.abspath(os.path.expanduser(allc_paths))
        if os.path.exists(input_path + '.tbi'):
            allc_paths = [allc_paths]
        else:
            df = pd.read_csv(input_path, sep='\t', header=None, usecols=[0])
            allc_paths = df.iloc[:, 0].tolist()
    allc_paths = [os.path.abspath(os.path.expanduser(allc_path))
                  for allc_path in allc_paths]
    outdir = os.path.abspath(os.path.expanduser(outdir))
    if not reference is None:
        reference = os.path.abspath(os.path.expanduser(reference))
        message = os.path.basename(reference)
        formats, columns, dimensions = ['H', 'H'], ['mc', 'cov'], ['chrom']
        usecols = [4, 5]
    else:
        message = ''
        formats, columns, dimensions = ['Q', 'H', 'H'], ['pos', 'mc', 'cov'], ['chrom']
        usecols = [1, 4, 5]
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    if not Path_to_chrom is None:
        Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
        df = pd.read_csv(Path_to_chrom, sep='\t', header=None, usecols=[0])
        chroms = df.iloc[:, 0].tolist()
    else:
        raise ValueError("Please provide chrom size path")
    manager = multiprocessing.Manager()
    queue1 = manager.Queue()
    pool = multiprocessing.Pool(n_jobs)
    watcher = pool.apply_async(mzs_merger, (formats, columns,
                                            dimensions, message, queue1))
    jobs = []
    for allc_path in allc_paths:
        outfile = os.path.join(outdir, os.path.basename(allc_path).replace(ext, '.mz'))
        if os.path.exists(outfile):
            print(f"{outfile} existed, skip.")
            continue
        outdir1 = outfile + '.tmp'
        if not os.path.exists(outdir1):
            os.mkdir(outdir1)
        tbi = pysam.TabixFile(allc_path)
        all_chroms = tbi.contigs
        tbi.close()
        # random.shuffle(chroms)
        chroms = [c for c in chroms if c in all_chroms]
        for chrom in chroms:
            if not reference is None:
                job = pool.apply_async(allc2mz_worker_with_ref,
                                       (allc_path, outdir1, reference, chrom, chroms,
                                        formats, columns, dimensions, usecols,
                                        missing_value, chunksize, pr, pa,
                                        sep, queue1))
            else:
                job = pool.apply_async(allc2mz_worker_without_ref,
                                       (allc_path, outdir1, chrom, chroms,
                                        formats, columns, dimensions,
                                        usecols, chunksize, sep, queue1))
            jobs.append(job)
    for job in jobs:
        r = job.get()
    queue1.put('kill')
    pool.close()
    pool.join()
    manager._process.terminate()
    manager.shutdown()


# ==========================================================
def _isCG(context):
    return context[:2] == b'CG'


# ==========================================================
def _isCH(context):
    return not _isCG(context)


# ==========================================================
def generate_context_ssi(Input, output=None, pattern='CGN'):
    if pattern == 'CGN':
        judge_func = _isCG
    else:  # CH
        judge_func = _isCH
    if output is None:
        output = Input + '.' + pattern + '.bmi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(Input)
    reader.category_ssi(output=output, formats=['I'], columns=['ID'],
                        dimensions=['chrom'], col=2, match_func=judge_func,
                        chunksize=2000)
    reader.close()


# ==========================================================
if __name__ == "__main__":
    pass
