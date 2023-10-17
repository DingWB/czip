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
from .bmz import Writer
# ==========================================================
class AllC:
    def __init__(self, Genome=None, Output="hg38_allc.mz",
                 pattern="C", jobs=4):
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
# =============================================================================
def print_allc(genome=None):
    """
    Print positions for all c in the reference genome. Example:
        python ~/Scripts/python/tbmate.py print_allc -g ~/genome/hg38/hg38.fa > hg38_allc.bed
    Parameters
    ----------
    genome: path
        reference genome.
    """
    genome = os.path.abspath(os.path.expanduser(genome))
    records = SeqIO.parse(genome, "fasta")
    ID = 0
    for record in records:
        chrom = record.id
        for i in range(record.seq.__len__() - 1):  # 0-based
            base = record.seq[i:i + 1].upper()
            if base.__str__() == 'C':  # forward strand
                context = record.seq[i: i + 3].upper().__str__()  # pos, left l1 base pair and right l2 base pair
                strand = '+'
            elif base.reverse_complement().__str__() == 'C':  # reverse strand
                context = record.seq[i - 2:i + 1].reverse_complement().upper().__str__()
                strand = '-'
            else:
                continue
            try:
                sys.stdout.write(f"{chrom}\t{i}\t{i + 1}\t{context}\t{ID}\t{strand}\n")
            except:
                sys.stdout.close()
                sys.exit()
            # position is 0-based (start) 1-based (end position, i+1)
            ID += 1
# =============================================================================
def view(Input="test_bed.bmzip"):
    reader = BmzReader(Input, 'rb')
    while reader._block_raw_length > 0:
        for r in struct.iter_unpack(reader.format, reader._buffer):
            yield r
        reader._load_block()
    reader.close()
# =============================================================================
if __name__ == "__main__":
    pass