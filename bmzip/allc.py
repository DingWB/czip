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
from multiprocessing import Pool
from Bio import SeqIO
# import pyximport
# pyximport.install(pyimport=True) #pyximport.install(pyimport=True)
from .utils import WriteC
from .bmz import Reader, Writer, allc2mz_mp
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
        pool = Pool(self.n_jobs)
        jobs = []
        for record in self.records:
            job = pool.apply_async(self.func, (record,self.outdir))
            jobs.append(job)
        for job in jobs:
            job.get()
        pool.close()
        pool.join()

    def merge(self):
        writer=Writer(Output=self.Output, Formats=['Q','c','3s'],
                      Columns=['pos', 'strand', 'context'],
                      Dimensions=['chrom'],message=self.genome)
        writer.catmz(Input=f"{self.outdir}/*.mz")

    def run(self):
        self.writePattern()
        self.merge()
        os.system(f"rm -rf {self.outdir}")
# ==========================================================
def allc2mz(allc_path, output, reference=None, n_jobs=12,
            chunksize=5000, verbose=0, path_to_chrom=None):
    """
    convert allc.tsv.gz to .mz file.

    Parameters
    ----------
    allc_path : path
        path to allc.tsv.gz, should has .tbi index.
    output : path
        path to .mz file
    reference : path
        path to reference coordinates.
    n_jobs : int
        number of cpu to process each allc file
    chunksize : int
        default is 5000
    verbose : int
    path_to_chrom : path
        path to chrom_size path or similar file containing chromosomes order,
        the first columns should be chromosomes, tab separated and no header.

    Returns
    -------

    """
    allc_path = os.path.abspath(os.path.expanduser(allc_path))
    if not os.path.exists(allc_path + '.tbi'):
        raise ValueError(f"allc file {allc_path} not found, please create .tbi index.")
    output = os.path.abspath(os.path.expanduser(output))
    outdir = output + '.tmp'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    if not reference is None:
        reference = os.path.abspath(os.path.expanduser(reference))
        formats, columns, dimensions = ['H', 'H'], ['mc', 'cov'], ['chrom']
        message = os.path.basename(reference)
        usecols = [4, 5]
    else:
        formats, columns, dimensions = ['Q', 'H', 'H'], ['pos', 'mc', 'cov'], ['chrom']
        message = ""
        usecols = [1, 4, 5]
    allc2mz_mp(allc_path, output, reference, message,
               formats, columns, dimensions, usecols, [0, 0],
               chunksize, 0, 1, n_jobs, '\t', verbose, None, path_to_chrom)


# ==========================================================
def __allcs2mzs_worker_ref(allc_path, outfile, path_to_chrom=None,
                           sep='\t', cols=[1, 4, 5], dtypes=[int, int, int],
                           columns=['pos', 'mc', 'cov'],
                           chunksize=100000):
    """
    convert allc.tsv.gz to .mz file.

    Parameters
    ----------
    allc_path : path
        path to allc.tsv.gz, should has .tbi index.
    output : path
        path to .mz file
    reference : path
        path to reference coordinates.
    n_jobs : int
        number of cpu to process each allc file
    chunksize : int
        default is 5000
    verbose : int
    path_to_chrom : path
        path to chrom_size path or similar file containing chromosomes order,
        the first columns should be chromosomes, tab separated and no header.

    Returns
    -------

    """
    global PosDict
    print(allc_path)

    def parse_record(line, cols):
        values = line.rstrip('\n').split(sep)
        return [values[i] for i in cols]

    if not os.path.exists(allc_path + '.tbi'):
        raise ValueError(f"allc file {allc_path} not found, please create .tbi index.")
    writer = Writer(outfile, Formats=['H', 'H'],
                    Columns=columns[1:], Dimensions=['chrom'])
    tbi = pysam.TabixFile(allc_path)
    chroms = tbi.contigs
    if not path_to_chrom is None:
        path_to_chrom = os.path.abspath(os.path.expanduser(path_to_chrom))
        df = pd.read_csv(path_to_chrom, sep='\t', header=None, usecols=[0])
        chroms = df.iloc[:, 0].tolist()
    for chrom in chroms:
        print(chrom, '\t' * 5, end='\r')
        df = pd.DataFrame([
            parse_record(line, cols) for line in tbi.fetch(chrom)],
            columns=columns
        )
        for col, dt in zip(columns, dtypes):
            df[col] = df[col].map(dt)
        df.set_index(columns[0], inplace=True)
        df = df.reindex(index=PosDict[tuple([chrom])], fill_value=0)
        while df.shape[0] > 0:
            data = df.iloc[:chunksize].apply(
                lambda x: struct.pack(f"<{writer.fmts}", *x.tolist()),
                axis=1).sum()
            writer.write_chunk(data, [chrom])
            try:
                df = df.iloc[chunksize:]
            except:
                break
    writer.close()
    tbi.close()


# ==========================================================
def load_reference(reference, dim=None, col=0):
    reader = Reader(reference)
    r = [record[col] for record in reader.fetch(dim)]
    reader.close()
    return dim, r


def allcs2mzs(allc_paths, outdir, reference=None, n_jobs=12,
              path_to_chrom=None, ext='.allc.tsv.gz'):
    if reference is None:
        raise ValueError("Converint multiple allc files without reference, please "
                         "use allc2mz and use n_jobs to do the parallely running.")
    outdir = os.path.abspath(os.path.expanduser(outdir))
    if isinstance(allc_paths, str):
        input_path = os.path.abspath(os.path.expanduser(allc_paths))
        df = pd.read_csv(input_path, sep='\t', header=None, usecols=[0])
        allc_paths = df.iloc[:, 0].tolist()
    allc_paths = [os.path.abspath(os.path.expanduser(allc_path))
                  for allc_path in allc_paths]
    print("Loading reference coordinates")
    refer_pos_dict = {}
    if n_jobs > 1:
        pool = Pool(n_jobs)
        jobs = []
        reader = Reader(reference)
        all_dims = reader.chunk_info.chunk_dims.tolist()
        reader.close()
        for dim in all_dims:
            job = pool.apply_async(load_reference, (reference, dim,))
            jobs.append(job)
        for job in jobs:
            r = job.get()
            refer_pos_dict[r[0]] = r[1]
        pool.close()
        pool.join()
    else:
        reader = Reader(reference)
        all_dims = reader.chunk_info.chunk_dims.tolist()
        for dim in all_dims:
            print(dim, end='\r')
            r = load_reference(reference, dim)
            refer_pos_dict[r[0]] = r[1]
        reader.close()

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    def init_worker(refer_pos_dict):
        global PosDict
        PosDict = refer_pos_dict

    pool = Pool(n_jobs, initializer=init_worker, initargs=(refer_pos_dict,))
    jobs = []
    for allc_path in allc_paths:
        outfile = os.path.join(outdir, os.path.basename(allc_path).replace(ext, '.mz'))
        job = pool.apply_async(__allcs2mzs_worker_ref,
                               (allc_path, outfile, path_to_chrom,))
        jobs.append(job)
    for job in jobs:
        r = job.get()
    pool.close()
    pool.join()


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
