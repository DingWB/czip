#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:34:38 2020

@author: DingWB
"""
import sys
import os
import struct
import pandas as pd
import pysam
import multiprocessing
from Bio import SeqIO
from collections import defaultdict
import time
# import pyximport
# pyximport.install(pyimport=True) #pyximport.install(pyimport=True)
from .utils import WriteC
from .bmz import Reader, Writer, get_dtfuncs
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
def allc2mz_worker_with_ref(allc_path, outdir, reference, chrom,
                            formats=['H', 'H'], columns=['mc', 'cov'],
                            dimensions=['chrom'], usecols=[4, 5],
                            missing_value=[0, 0], chunksize=20000,
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
    ref_records = ref_reader.fetch(tuple([chrom]))
    data = b''
    na_value_bytes = struct.pack(f"<{writer.fmts}", *missing_value)
    i = 0
    while True:
        try:
            row_query = records.__next__().rstrip('\n').split(sep)
        except:
            break
        row_ref = ref_records.__next__()
        while row_ref[pr] < int(row_query[pa]):
            data += na_value_bytes
            i += 1
            try:
                row_ref = ref_records.__next__()
            except:
                break
        if row_ref[pr] == int(row_query[pa]):  # match
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
        q.put(tuple([outdir, chrom]))
    return outdir, chrom


# ==========================================================
def allc2mz_worker_without_ref(allc_path, outdir, chrom,
                               formats=['Q', 'H', 'H'], columns=['pos', 'mc', 'cov'],
                               dimensions=['chrom'], usecols=[1, 4, 5],
                               chunksize=5000, sep='\t', q=None):
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
        q.put(tuple([outdir, chrom]))
    return outdir, chrom


def mzs_merger(chroms, formats, columns, dimensions, message, q):
    N = len(chroms)
    finished_jobs = {}
    while 1:
        result = q.get()  # if q is empty, it will wait.
        if result == 'kill':
            print("Done!")
            break
        outdir, chrom = result
        if outdir not in finished_jobs:
            finished_jobs[outdir] = [chrom]
        else:
            finished_jobs[outdir].append(chrom)
        keys = list(finished_jobs.keys())
        for outdir in keys:
            if len(finished_jobs[outdir]) < N:
                continue  # wait
            outfile = outdir[:-4]  # .tmp
            writer = Writer(Output=outfile, Formats=formats,
                            Columns=columns, Dimensions=dimensions,
                            message=message)
            writer.catmz(Input=f"{outdir}/*")
            os.system(f"rm -rf {outdir}")
            basename = os.path.basename(outfile)
            print(f"{basename} finished", "\t" * 4, end='\t')
            del finished_jobs[outdir]  # finished one allc, remove from monitoring
    return


# ==========================================================
def allc2mz(allc_paths, outdir, reference, missing_value=[0, 0],
            chunksize=5000, pr=0, pa=1,
            n_jobs=8, sep='\t', Path_to_chrom=None,
            ext='.allc.tsv.gz'):
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
    watcher = pool.apply_async(mzs_merger, (chroms, formats, columns,
                                            dimensions, message, queue1))
    jobs = []
    for allc_path in allc_paths:
        outfile = os.path.join(outdir, os.path.basename(allc_path).replace(ext, '.mz'))
        outdir1 = outfile + '.tmp'
        if not os.path.exists(outdir1):
            os.mkdir(outdir1)
        for chrom in chroms:
            if not reference is None:
                job = pool.apply_async(allc2mz_worker_with_ref,
                                       (allc_path, outdir1, reference, chrom,
                                        formats, columns, dimensions, usecols,
                                        missing_value, chunksize, pr, pa,
                                        sep, queue1))
            else:
                job = pool.apply_async(allc2mz_worker_without_ref,
                                       (allc_path, outdir1, chrom,
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
    return context[:2] == 'CG'


# ==========================================================
def _isCH(context):
    return not _isCG(context)


# ==========================================================
def generate_context_ssi(Input, output=None, formats=['I'], columns=['ID'],
                         dimensions=['chrom'], col=2, pattern='CGN',
                         chunksize=5000):
    if pattern == 'CGN':
        judge_func = _isCG
    else:  # CH
        judge_func = _isCH
    if output is None:
        output = Input + '.' + pattern + '.bmi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(Input)
    writer = Writer(output, Formats=formats, Columns=columns,
                    Dimensions=dimensions, fileobj=None,
                    message=Input)
    data = b''
    for dim in reader.dim2chunk_start:
        for i, record in enumerate(reader.fetch(dim)):
            if judge_func(record[col]):
                data += struct.pack(f"<{writer.fmts}", i + 1)
            if ((i + 1) % chunksize) == 0 and len(data) > 0:
                writer.write_chunk(data, dim)
                data = b''
        if len(data) > 0:
            writer.write_chunk(data, dim)
            data = b''
    writer.close()


# ==========================================================
if __name__ == "__main__":
    pass
