#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 17:34:38 2020

@author: DingWB
"""
import sys
import os
import struct
import gzip
import pandas as pd
from multiprocessing import Pool
from Bio import SeqIO
import pyximport
pyximport.install(pyimport=True) #pyximport.install(pyimport=True)
from .utils import WriteC
from .bmz import Writer
# =============================================================================
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
        self.chroms=[]
        for record in self.records:
            job = pool.apply_async(self.func, (record,self.outdir))
            self.chroms.append(record.id)
            jobs.append(job)
        # self.results=[]
        for job in jobs:
            # self.results.append(job.get())
            job.get()
        pool.close()
        pool.join()

    def merge(self):
        writer=Writer(Output=self.Output, Formats=['Q','c','3s'],
                      Names=['pos', 'strand', 'context'],
                      Tags=['chrom'],message=self.genome)
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

def open1(infile):
    if infile.endswith('.gz'):
        f=gzip.open(infile,'rb')
    else:
        f=open(infile,'rb')
    return f
# =============================================================================
# class viewer: #python ~/Scripts/python/tbmate.py viewer -- --interactive
# 	def __init__(self, x=None,input=None,columns=['mc','cov']): #vw=viewer(x="~/genome/hg38/annotations/hg38_allc.bed.gz",input="test*.bin.gz")
# 		"""
# 		View allc.bin.gz, for example:
# 			tbmate viewer -x ~/genome/hg38/annotations/hg38_allc.bed.gz -i "test*.bin.gz"  view |head
# 		Parameters
# 		----------
# 		x: path
# 			path for the index
# 		input: path
# 			path of input *.bin or *.bin.gz files, it should be quoted if wildcard is included,
# 			for example: --input "test/*.bin.gz"
# 		"""
# 		self.x=os.path.abspath(os.path.expanduser(x))
# 		self.input=glob.glob(os.path.abspath(os.path.expanduser(input)))
# 		self.columns=columns
# 		self.get_header()
# 		self.prepare_input()
#
# 	def get_header(self):
# 		self.header_dict={}
# 		for infile in self.input:
# 			sname=os.path.basename(infile).split('.')[0]
# 			self.header_dict[sname]=header(infile,ret=True) #N, dtypes, idx
#
# 	def prepare_input(self):
# 		# print("Preparing the input files..")
# 		self.fi={}
# 		for infile in self.input:
# 			sname=os.path.basename(infile).split('.')[0]
# 			self.fi[sname]=open1(infile)
#
# 	def header(self):
# 		sys.stdout.write(f"sampleID\tN\tdtypes\tidx\n")
# 		for sname in self.header_dict:
# 			sys.stdout.write(sname)
# 			for v in self.header_dict[sname]:
# 				sys.stdout.write(f"\t{v}")
# 			sys.stdout.write('\n')
#
# 	def read(self,sname,chunksize):
# 		self.fi[sname].seek(HEADER_LEN)
# 		dtypes = self.header_dict[sname][1]
# 		fmts = ''.join([dtype_fmt[dtype] for dtype in dtypes])
# 		sizes = [struct.calcsize(fmt) for fmt in fmts]
# 		s = sum(sizes)
# 		columns=[f"{sname}."+str(col) for col in self.columns]
# 		read_size=chunksize * s
# 		while True:
# 			R=[]
# 			for i in struct.iter_unpack(fmts,self.fi[sname].read(read_size)):
# 				R.append(i)
# 			df1=pd.DataFrame(R,columns=columns)
# 			yield df1
#
# 	def get_chunk(self, chunksize):
# 		data={}
# 		for infile in self.input:
# 			sname = os.path.basename(infile).split('.')[0]
# 			data[sname]=self.read(sname,chunksize)
# 		sys.stdout.write('\t'.join(['chrom', 'pos', 'context', 'strand']))
# 		for sname in data:
# 			for col in self.columns:
# 				sys.stdout.write("\t"+str(sname)+'.'+str(col))
# 		sys.stdout.write("\n")
# 		for df in pd.read_csv(self.x, sep='\t', usecols=[0, 2, 3, 5], header=None,
# 							  chunksize=chunksize, names=['chrom', 'pos', 'contet', 'strand']):
# 			for sname in data:
# 				df1=next(data[sname])
# 				for col in df1:
# 					df[col]=df1[col]
# 			yield df
#
# 	def close(self):
# 		for sname in self.fi:
# 			self.fi[sname].close()
#
# 	def view(self,chunksize=100000):
# 		for df in self.get_chunk(chunksize=chunksize):
# 			try:
# 				for values in df.applymap(str).values:
# 					sys.stdout.write('\t'.join(values)+'\n')
# 			except:
# 				self.close()
# 				sys.stdout.close()
# 				sys.exit()
# 		self.close()
#
# 	def sread(self,sname,start,L):
# 		dtypes = self.header_dict[sname][1]
# 		fmts = ''.join([dtype_fmt[dtype] for dtype in dtypes])
# 		sizes = [struct.calcsize(fmt) for fmt in fmts]
# 		s = sum(sizes)
# 		read_size = L * s
# 		self.fi[sname].seek(start*s+HEADER_LEN)
# 		columns=[f"{sname}."+str(col) for col in self.columns]
# 		R=[]
# 		for i in struct.iter_unpack(fmts,self.fi[sname].read(read_size)):
# 			R.append(i)
# 		df1=pd.DataFrame(R,columns=columns)
# 		return df1
#
# 	def query(self,regions): # regions="chr1:72135220-72135258"
# 		if type(regions)==str:
# 			regions=[regions]
# 		sys.stdout.write('\t'.join(['chrom', 'pos', 'context', 'strand']))
# 		for sname in self.fi:
# 			for col in self.columns:
# 				sys.stdout.write("\t" + str(sname) + '.' + str(col))
# 		sys.stdout.write("\n")
# 		tb = tabix.open(self.x)
# 		for region in sorted(regions):
# 			records=tb.querys(region) #query(seqname,start,end)
# 			R=[]
# 			for record in records:
# 				R.append(record)
# 			df=pd.DataFrame(R,columns=['chrom','begin','end','context','ID','strand'])
# 			start=int(df.ID.iloc[0])
# 			L=df.shape[0]
# 			df=df.loc[:,['chrom','end','context','strand']]
# 			for sname in self.fi:
# 				df1 = self.sread(sname,start,L)
# 				for col in df1:
# 					df[col] = df1[col].map(str)
# 				# print(df.shape,df.dtypes)
# 			try:
# 				for values in df.values:
# 					sys.stdout.write('\t'.join(values) + '\n')
# 			except:
# 				self.close()
# 				sys.stdout.close()
# 				sys.exit()
# 		self.close()
# =============================================================================
def create_index(infile, skipn=0, sep='\t',save=False):
    f = open1(infile)
    for i in range(skipn):
        line = f.readline()
    line = f.readline()
    if isinstance(line, bytes):
        line = line.decode('utf-8')
    line_num = 0
    chrom = ""
    start = 0
    R = []
    r = []
    while line:
        chr = line.strip().split(sep)[0]
        if chr != chrom:
            start = line_num
            if len(r) == 0:
                r.extend([chr, start])
            elif len(r) == 2:
                r.append(start - 1)
                R.append(r)
                r = [chr, start]
            chrom = chr
            # if len(R)>1:
            #     print(chrom,r,R[-1])
            # else:
            #     print(chrom,r)
        line = f.readline()
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        line_num += 1
    f.close()
    if save:
        df = pd.DataFrame(R)
        df.to_csv(infile + '.fai', sep='\t', index=False, header=False)
    return R
# =============================================================================
def get_chrom_dict(infile=None,sep='\t',skipn=0,save=True):
    if not os.path.exists(os.path.expanduser(infile+'.fai')):
        R=create_index(infile,skipn=skipn,sep=sep,save=save)
    else:
        df=pd.read_csv(infile+'.fai',sep='\t',header=None)
        R=df.tolist()
    chrom_dict={}
    for r in R:
        chrom_dict[r[0]]=tuple(r[1:])
    return chrom_dict
# =============================================================================
def read1(f,start,fmt):
    f.seek(start)
    r = f.read(struct.calcsize(fmt))
    return struct.unpack(fmt,r)[0]


if __name__ == "__main__":
    pass