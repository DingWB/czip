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
from multiprocessing import Pool
from Bio import SeqIO
# import pyximport
# pyximport.install(pyimport=True) #pyximport.install(pyimport=True)
from .utils import WriteC
from .bmz import Writer, allc2mz_mp
# ==========================================================
class AllC:
    def __init__(self, Genome=None, Output="hg38_allc.mz",
                 pattern="C", jobs=12):
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
        self.n_jobs=jobs if not jobs is None else os.cpu_count()
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
def allc2mz(allc_path, output, reference=None, jobs=12,
            chunksize=5000, verbose=0):
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
               chunksize, 0, 1, jobs, '\t', verbose)
# ==========================================================
if __name__ == "__main__":
    pass