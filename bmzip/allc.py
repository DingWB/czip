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
import random
# from .utils import WriteC
# sys.path.append(os.path.dirname(__file__))
from .bmz import Reader, Writer, get_dtfuncs


# ==========================================================
def WriteC(record, outdir, chunksize=5000):
    # cdef int i, N
    # cdef char* chrom, base, context, strand
    chrom = record.id
    outfile = os.path.join(outdir, chrom + ".mz")
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
    def __init__(self, Genome=None, Output="hg38_allc.mz",
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
        writer.catmz(Input=f"{self.outdir}/*.mz")

    def run(self):
        self.writePattern()
        self.merge()
        if not self.keep_temp:
            os.system(f"rm -rf {self.outdir}")

def allc2mz(allc_path, outfile, reference=None, missing_value=[0, 0],
            pr=0, pa=1, sep='\t', Path_to_chrom=None, chunksize=2000):
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
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    allc_path = os.path.abspath(os.path.expanduser(allc_path))
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
        formats, columns, dimensions = ['H', 'H'], ['mc', 'cov'], ['chrom']
        usecols = [4, 5]
    else:
        message = ''
        formats, columns, dimensions = ['Q', 'H', 'H'], ['pos', 'mc', 'cov'], ['chrom']
        usecols = [1, 4, 5]
    writer = Writer(outfile, Formats=formats, Columns=columns,
                    Dimensions=dimensions, message=message)
    dtfuncs = get_dtfuncs(formats, tobytes=False)

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
def generate_ssi(Input, output=None, pattern="+CGN"):
    """
    Generate ssi (subset index) for a given input .mz

    Parameters
    ----------
    Input : .mz
    output : .bmi
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
        output = Input + '.' + pattern + '.bmi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(Input)
    reader.category_ssi(output=output, formats=['I'], columns=['ID'],
                        dimensions=['chrom'], match_func=judge_func,
                        chunksize=2000)
    reader.close()

def prepare_sky(smk=None, sky=None, indir=None, outdir=None, allc_path=None,
                reference=None, ref_prefix=None, chrom=None, chrom_prefix=None,
                gcp=False, bucket=None, cpu=24, name=None):
    if smk is None:
        smk = os.path.join(os.path.dirname(__file__),
                           "data/snakemake_template/run_allc2mz.smk")
    if sky is None:
        sky = os.path.join(os.path.dirname(__file__),
                           "data/skypilot_template/run_allc2mz.yaml")
    if name is None:
        name = 'allc2mz'
    workdir = os.path.basename("./")
    D = {
        'indir': indir, 'outdir': outdir, 'allc_path': allc_path,
        'reference': reference, 'ref_prefix': ref_prefix, 'chrom': chrom,
        'chrom_prefix': chrom_prefix, 'gcp': gcp, 'bucket': bucket,
        'cpu': cpu
    }
    for k in D:
        if D[k] is None:
            D[k] = ''
        else:
            if k not in ['bucket', 'cpu']:
                value = D[k]
                D[k] = f"{k}={value}"
    # D['smk']=smk
    D['name'] = name
    D['workdir'] = workdir
    with open(sky, 'r') as f:
        template = f.read()
    # with open(out_yaml,'w') as f:
    #     f.write(template.format(**D))
    print(template.format(**D))
    print("# sky launch -c test 1.yaml")
    print("# sky spot launch -y -n job job.yaml")


def mzs_merger(outfile, formats, columns, dimensions, message, n_batch, q):
    finished_jobs = {}
    N = {}
    chunk_id_tmp = 0
    writer = Writer(Output=outfile, Formats=formats,
                    Columns=columns, Dimensions=dimensions,
                    message=message)
    dtfuncs = get_dtfuncs(writer.Formats)
    while 1:
        result = q.get()  # if q is empty, it will wait.
        if result == 'kill':
            print("Done!")
            writer.close()
            break
        dim, chunk_id, data = result  # chunk_id==''last_one''
        if dim not in finished_jobs:
            print(dim)
            finished_jobs[dim] = {}
            N[dim] = {}
        if chunk_id not in finished_jobs[dim]:
            finished_jobs[dim][chunk_id] = data
            N[dim][chunk_id] = 0
            # print(dim,chunk_id,end='\r')
        else:
            finished_jobs[dim][chunk_id] += data  # sum this chunk for another 100 files
            N[dim][chunk_id] += 1  # record how many batch have been processed
        if N[dim][chunk_id_tmp] == n_batch:
            data = b''.join([struct.pack(f"<{writer.fmts}",
                                         *[func(v) for v, func in zip(values, dtfuncs)])
                             for values in finished_jobs[dim][chunk_id]])
            print(dim, chunk_id_tmp, N[dim][chunk_id_tmp], n_batch)
            writer.write_chunk(data, dim)
            del finished_jobs[dim][chunk_id]  # finished one allc, remove from monitoring
            chunk_id_tmp += 1
            if chunk_id_tmp not in N[dim]:
                chunk_id_tmp = 0


def merge_mz_worker(outfile_cat, outdir, chrom, dims, formats, chunksize=5000):
    # print(chrom, "started",end='\r')
    outname = os.path.join(outdir, chrom + '.mz')
    reader = Reader(outfile_cat)
    data = None
    for dim in dims:  # each dim is a file of the same chrom
        values = np.array([record for record in reader.__fetch__(dims=dim, s=0, e=2)])
        if data is None:
            data = values.copy()
        else:
            data += values
        # q.put(tuple([dim, chunk_id, data]))
    reader.close()
    writer = Writer(outname, Formats=formats,
                    Columns=reader.header['Columns'],
                    Dimensions=reader.header['Dimensions'][:1],
                    message='outfile_cat')
    byte_data, i = b'', 0
    for values in data:
        byte_data += struct.pack(f"<{writer.fmts}",
                                 *values)
        i += 1
        if i > chunksize:
            writer.write_chunk(byte_data, [chrom])
            byte_data, i = b'', 0
    if len(byte_data) > 0:
        writer.write_chunk(byte_data, [chrom])
    writer.close()
    print(chrom, "done")
    return


def merge_mz(indir="/anvil/scratch/x-wding2/Projects/mouse-pfc/data/pseudo_cell/mz-CGN",
             mz_paths=None, outfile="merged.mz", n_jobs=12, formats=['I', 'I'],
             Path_to_chrom="~/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt",
             keep_cat=False):
    outfile = os.path.abspath(os.path.expanduser(outfile))
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    if mz_paths is None:
        mz_paths = os.listdir(indir)
    reader = Reader(os.path.join(indir, mz_paths[0]))
    header = reader.header
    reader.close()
    outfile_cat = outfile + '.cat.mz'
    writer = Writer(Output=outfile_cat, Formats=header['Formats'],
                    Columns=header['Columns'], Dimensions=header['Dimensions'],
                    message="catmz")
    writer.catmz(Input=[os.path.join(indir, mz_path) for mz_path in mz_paths],
                 add_dim=True)
    Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
    df = pd.read_csv(Path_to_chrom, sep='\t', header=None, usecols=[0])
    chroms = df.iloc[:, 0].tolist()

    reader = Reader(outfile_cat)
    chrom_col = reader.header['Dimensions'][0]
    chunk_info = reader.chunk_info
    reader.close()
    # manager = multiprocessing.Manager()
    # queue1 = manager.Queue()
    pool = multiprocessing.Pool(n_jobs)
    # watcher = pool.apply_async(mzs_merger, (outfile, formats, columns,
    #                                         dimensions, message, n_batch, queue1))
    jobs = []
    outdir = outfile + '.tmp'
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    for chrom in chroms:
        dims = chunk_info.loc[chunk_info[chrom_col] == chrom].index.tolist()
        if len(dims) == 0:
            continue
        job = pool.apply_async(merge_mz_worker,
                               (outfile_cat, outdir, chrom, dims, formats))
        jobs.append(job)
    for job in jobs:
        r = job.get()
    # queue1.put('kill')
    pool.close()
    pool.join()
    # manager._process.terminate()
    # manager.shutdown()
    if not keep_cat:
        os.remove(outfile_cat)
    writer = Writer(Output=outfile, Formats=formats,
                    Columns=header['Columns'], Dimensions=header['Dimensions'],
                    message="catmz")
    writer.catmz(Input=f"{outdir}/*.mz")
    os.system(f"rm -rf {outdir}")


# ==========================================================
def extractCG(Input=None, outfile=None, bmi=None, chunksize=5000,
              merge_strand=True):
    """

    Parameters
    ----------
    mz_path :
    outfile :
    bmi : path
        bmi should be bmi to mm10_with_chrL.allc.mz.CGN.bmi, not forward
        strand bmi, but after merge (if merge_strand is True), forward bmi
        mm10_with_chrL.allc.mz.+CGN.bmi should be used to generate
         reference, one can
        run: bmzip extract -m mm10_with_chrL.allc.mz
        -o mm10_with_chrL.allCG.forward.mz
        -b mm10_with_chrL.allc.mz.+CGN.bmi and use
        mm10_with_chrL.allCG.forward.mz as new reference.
    chunksize :int
    merge_strand: bool
        after merging, only forward strand would be kept, reverse strand values
        would be added to the corresponding forward strand.

    Returns
    -------

    """
    mz_path = os.path.abspath(os.path.expanduser(Input))
    bmi_path = os.path.abspath(os.path.expanduser(bmi))
    bmi_reader = Reader(bmi_path)
    reader = Reader(mz_path)
    writer = Writer(outfile, Formats=reader.header['Formats'],
                    Columns=reader.header['Columns'],
                    Dimensions=reader.header['Dimensions'],
                    message=bmi_path)
    dtfuncs = get_dtfuncs(writer.Formats)
    for dim in reader.dim2chunk_start.keys():
        print(dim)
        IDs = bmi_reader.get_ids_from_bmi(dim)
        if len(IDs.shape) != 1:
            raise ValueError("Only support 1D bmi now!")
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
    bmi_reader.close()


# ==========================================================
if __name__ == "__main__":
    import fire

    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire()
