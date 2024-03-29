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
from scipy import stats
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
        L = len(context)
        if L < 3:
            if L == 0:
                context = "CNN"
            else:
                context = context + 'N' * (3 - L)

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
           Formats=['B', 'B'], Columns=['mc', 'cov'], Dimensions=['chrom'],
           usecols=[4, 5], pr=0, pa=1, sep='\t', Path_to_chrom=None,
           chunksize=5000):
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
def generate_ssi1(input, output=None, pattern="CGN"):
    """
    Generate 1D ssi (subset index) for a given input .cz, 1D means calculating
    the ID list for a given pattern.

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
        output = input + '.' + pattern + '.ssi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(input)
    reader.category_ssi(output=output, formats=['I'], columns=['ID'],
                        dimensions=['chrom'], match_func=judge_func,
                        chunksize=2000)
    reader.close()


def generate_ssi2(input, output=None, bed=None,
                  n_jobs=4):  # 2D ssi
    """
    Generate subset index for a genomic region bed file. For example:
        czip generate_ssi2 -i ~/Ref/mm10/annotations/mm10_with_chrL.allc.cz \
        -o mm10_with_chrL.allc.genes_flank2k.ssi -b genes_flank2k.bed.gz -n 4
    Parameters
    ----------
    input :
    output :
    bed :
    n_jobs :

    Returns
    -------

    """
    bed = os.path.abspath(os.path.expanduser(bed))
    input = os.path.abspath(os.path.expanduser(input))
    if output is None:
        output = input + '.' + os.path.basename(bed) + '.ssi'
    else:
        output = os.path.abspath(os.path.expanduser(output))
    reader = Reader(input)
    reader.regions_ssi(output=output, bed=bed, n_jobs=n_jobs)
    reader.close()

# ==========================================================
def _fisher_worker(df):
    import warnings
    warnings.filterwarnings("ignore")
    from fast_fisher import fast_fisher_exact, odds_ratio
    # one verse rest
    columns = df.columns.tolist()
    snames = [col[:-3] for col in columns if col.endswith('.mc')]
    df['mc_sum'] = df.loc[:, [name for name in columns if name.endswith('.mc')]].sum(axis=1)
    df['cov_sum'] = df.loc[:, [name for name in columns if name.endswith('.cov')]].sum(axis=1)
    df['uc_sum'] = df.cov_sum - df.mc_sum

    def cal_fisher_or_p(x):
        uc = int(x[f"{sname}.cov"] - x[f"{sname}.mc"])
        a, b, c, d = int(x[f"{sname}.mc"]), uc, int(x.mc_sum - x[f"{sname}.mc"]), int(x.uc_sum - uc)
        Or = odds_ratio(a, b, c, d)
        Pval = fast_fisher_exact(a, b, c, d)
        return tuple(['%.3g' % Or, '%.3g' % Pval])

    for sname in snames:
        df[sname] = df.apply(cal_fisher_or_p, axis=1)
        df[f"{sname}.odd_ratio"] = df[sname].apply(lambda x: x[0])
        df[f"{sname}.pval"] = df[sname].apply(lambda x: x[1])
        df.drop([f"{sname}.cov", f"{sname}.mc", sname], axis=1, inplace=True)
    usecols = []
    for sname in snames:
        usecols.extend([f"{sname}.odd_ratio", f"{sname}.pval"])
    return df.reindex(columns=usecols)


# ==========================================================
def merge_cz_worker(outfile_cat, outdir, chrom, dims, formats,
                    block_idx_start, batch_nblock, chunksize=5000):
    # print(chrom, "started",end='\r')
    if formats in ['fraction', '2D', 'fisher']:
        ext = 'txt'
    else:
        ext = 'cz'
    outname = os.path.join(outdir, chrom + f'.{block_idx_start}.{ext}')
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
        if formats == 'fraction':
            values = np.array([[0] if v[1] == 0 else ['%.3g' % (v[0] / v[1])] for v in values])
        if data is None:
            data = values.copy()
        else:
            if formats in ["fraction", "2D", 'fisher']:
                data = np.hstack((data, values))
            else:
                data += values
        # q.put(tuple([dim, chunk_id, data]))
    if formats in ['fraction', '2D', 'fisher']:
        snames = [dim[1] for dim in dims]
        if formats == 'fraction':
            columns = snames
        else:
            columns = []
            for sname in snames:
                columns.extend([sname + '.mc', sname + '.cov'])
        df = pd.DataFrame(data, columns=columns)
        if formats == 'fisher':
            df = _fisher_worker(df)
        df.to_csv(outname, sep='\t', index=False)
        # print(chrom, block_idx_start, "done")
        return

    writer1 = Writer(outname, Formats=formats,
                     Columns=reader1.header['Columns'],
                     Dimensions=reader1.header['Dimensions'][:1],
                     message=outfile_cat)
    byte_data, i = b'', 0
    dtfuncs = get_dtfuncs(writer1.Formats)
    for values in data.tolist():
        byte_data += struct.pack(f"<{writer1.fmts}",
                                 *[func(v) for v, func in zip(values, dtfuncs)])
        i += 1
        if i > chunksize:
            writer1.write_chunk(byte_data, [chrom])
            byte_data, i = b'', 0
    if len(byte_data) > 0:
        writer1.write_chunk(byte_data, [chrom])
    writer1.close()
    reader1.close()
    # print(chrom, block_idx_start, "done", "\t" * 8, end='\r')
    return


def catchr(outdir, chrom, ext, batch_nblock, chunksize):
    outname = os.path.join(outdir, f"{chrom}.{ext}")
    block_idx_start = 0
    infile = os.path.join(outdir, chrom + f'.{block_idx_start}.{ext}')
    while os.path.exists(infile):
        for df in pd.read_csv(infile, sep='\t', chunksize=chunksize):
            if not os.path.exists(outname):
                df.to_csv(outname, sep='\t', index=False, header=True)
            else:
                df.to_csv(outname, sep='\t', index=False, header=False, mode='a')
        block_idx_start += batch_nblock
        infile = os.path.join(outdir, chrom + f'.{block_idx_start}.{ext}')
    return


def merge_cz(indir=None, cz_paths=None, class_table=None,
             outfile=None, prefix=None, n_jobs=12, formats=['H', 'H'],
             Path_to_chrom=None, reference=None,
             keep_cat=False, batchsize=10, temp=False, bgzip=True,
             chunksize=50000, ext='.cz'):
    """
    Merge multiple .cz files. For example:
    czip merge_cz -i ./ -o major_type.2D.txt -n 96 -f 2D \
                          -P ~/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt \
                          -r ~/Ref/mm10/annotations/mm10_with_chrL.allCG.forward.cz

    Parameters
    ----------
    indir :path
        If cz_paths is not provided, indir will be used to get cz_paths.
    cz_paths :paths
    class_table: path
        If class_table is given, multiple output will be generated based on the
        snames and class from this class_table, each output will have a suffix of
        class name in this table.
    outfile : path
    n_jobs :int
    formats : str of list
        Could be fraction, 2D, fisher or list of formats.
        if formats is a list, then mc and cov will be summed up and write to .cz file.
        otherwise, if formats=='fraction', summed mc divided by summed cov
        will be calculated and written to .txt file. If formats=='2D', mc and cov
        will be kept and write to .txt matrix file.
    Path_to_chrom : path
        path to chrom size file.
    reference : path
        path to reference .cz file, only need if fraction="fraction" or "2D".
    keep_cat : bool
    batchsize :int
    temp : bool
    bgzip : bool
    chunksize : int

    Returns
    -------

    """
    if not class_table is None:
        df_class = pd.read_csv(class_table, sep='\t', header=None,
                               names=['sname', 'cell_class'])
        snames = [file.replace(ext, '') for file in os.listdir(indir)]
        df_class = df_class.loc[df_class.sname.isin(snames)]
        D = df_class.groupby('cell_class').sname.apply(
            lambda x: x.tolist()).to_dict()
        for key in D:
            print(key)
            cz_paths = [sname + ext for sname in D[key]]
            merge_cz(indir, cz_paths, class_table=None,
                     outfile=None, prefix=f"{prefix}.{key}", n_jobs=n_jobs,
                     formats=formats, Path_to_chrom=Path_to_chrom,
                     reference=reference, keep_cat=keep_cat,
                     batchsize=batchsize, temp=temp, bgzip=bgzip,
                     chunksize=chunksize, ext=ext)
        return None
    if outfile is None:
        if prefix is None:
            outfile = 'merged.cz' if formats not in ['fraction', '2D', 'fisher'] else 'merged.txt'
        else:
            outfile = f'{prefix}.cz' if formats not in ['fraction', '2D', 'fisher'] else f'{prefix}.txt'
    print(outfile)
    outfile = os.path.abspath(os.path.expanduser(outfile))
    if os.path.exists(outfile):
        print(f"{outfile} existed, skip.")
        return
    if cz_paths is None:
        cz_paths = [file for file in os.listdir(indir) if file.endswith(ext)]
    reader = Reader(os.path.join(indir, cz_paths[0]))
    header = reader.header
    reader.close()
    outfile_cat = outfile + '.cat.cz'
    # cat all .cz files into one .cz file, add a dimension to chunk (filename)
    writer = Writer(Output=outfile_cat, Formats=header['Formats'],
                    Columns=header['Columns'], Dimensions=header['Dimensions'],
                    message="catcz")
    writer.catcz(Input=[os.path.join(indir, cz_path) for cz_path in cz_paths],
                 add_dim=True)

    reader = Reader(outfile_cat)
    chrom_col = reader.header['Dimensions'][0]
    chunk_info = reader.chunk_info
    reader.close()

    # get chromosomes order
    input_chroms = chunk_info[chrom_col].unique().tolist()
    if not Path_to_chrom is None:
        Path_to_chrom = os.path.abspath(os.path.expanduser(Path_to_chrom))
        df = pd.read_csv(Path_to_chrom, sep='\t', header=None, usecols=[0])
        chroms = [chrom for chrom in df.iloc[:, 0].tolist() if chrom in input_chroms]
    else:
        chroms = sorted(input_chroms)
    # chrom_dims_dict=chunk_info.reset_index().groupby('chrom'
    #                                                    ).chunk_dims.apply(
    #     lambda x:x.unique().tolist()).to_dict()
    chrom_nblocks = chunk_info.reset_index().loc[:, [chrom_col, 'chunk_nblocks']
                    ].drop_duplicates().set_index(chrom_col).chunk_nblocks.to_dict()
    # chunk_info.set_index(chrom_col)
    # how many blocks can be multiplied by self.unit_size
    unit_nblock = int(writer._unit_size / (math.gcd(writer._unit_size, _BLOCK_MAX_LEN)))
    nunit_perbatch = int(np.ceil((chunk_info.chunk_nblocks.max() / batchsize
                                  ) / unit_nblock))
    batch_nblock = nunit_perbatch * unit_nblock  # how many block for each batch
    pool = multiprocessing.Pool(n_jobs)
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
    pool.close()
    pool.join()

    # First, merge different batch for each chrom
    if formats in ['fraction', '2D', 'fisher']:
        out_ext = 'txt'
    else:
        out_ext = 'cz'
    if out_ext == 'cz':
        for chrom in chroms:
            # merge batches into chrom (chunk)
            outname = os.path.join(outdir, f"{chrom}.{out_ext}")
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
            infile = os.path.join(outdir, chrom + f'.{block_idx_start}.{out_ext}')
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
                reader.close()
                block_idx_start += batch_nblock
                infile = os.path.join(outdir, chrom + f'.{block_idx_start}.{out_ext}')
            # write chunk tail
            writer.close()
    else:  # txt
        pool = multiprocessing.Pool(n_jobs)
        jobs = []
        for chrom in chroms:
            job = pool.apply_async(catchr,
                                   (outdir, chrom, out_ext, batch_nblock, chunksize))
            jobs.append(job)
        for job in jobs:
            r = job.get()
        pool.close()
        pool.join()

    # Second, merge chromosomes to outfile
    if out_ext == 'cz':  # merge chroms into final output
        writer = Writer(Output=outfile, Formats=formats,
                        Columns=header['Columns'], Dimensions=header['Dimensions'],
                        message="merged")
        writer.catcz(Input=[f"{outdir}/{chrom}.cz" for chrom in chroms])
    else:  # txt
        filenames = chunk_info.filename.unique().tolist()
        if formats == 'fraction':
            columns = filenames
        elif formats == '2D':
            columns = []
            for sname in filenames:
                columns.extend([sname + '.mc', sname + '.cov'])
        else:  # fisher
            columns = []
            for sname in filenames:
                columns.extend([sname + '.odd_ratio', sname + '.pval'])
        if not reference is None:
            reference = os.path.abspath(os.path.expanduser(reference))
            reader = Reader(reference)
        print("Merging chromosomes..")
        for chrom in chroms:
            print(chrom, "\t" * 4, end='\r')
            infile = os.path.join(outdir, f"{chrom}.{out_ext}")
            if not reference is None:
                df_ref = pd.DataFrame([
                    record for record in reader.fetch(tuple([chrom]))
                ], columns=reader.header['Columns'])
                # insert a column 'start'
                df_ref.insert(0, chrom_col, chrom)
                df_ref.insert(1, 'start', df_ref.iloc[:, 1].map(int) - 1)
                usecols = df_ref.columns.tolist() + columns
            # df=pd.read_csv(infile,sep='\t')
            for df in pd.read_csv(infile, sep='\t', chunksize=chunksize):
                if not reference is None:
                    df = pd.concat([df_ref.iloc[:chunksize].reset_index(drop=True),
                                    df.reset_index(drop=True)], axis=1)
                    df_ref = df_ref.iloc[chunksize:]
                if not os.path.exists(outfile):
                    df.reindex(columns=usecols).to_csv(outfile, sep='\t', index=False, header=True)
                else:
                    df.reindex(columns=usecols).to_csv(outfile, sep='\t', index=False,
                                                       header=False, mode='a')
        if not reference is None:
            reader.close()
    if not keep_cat:
        os.remove(outfile_cat)
    if not temp:
        print(f"Removing temp dir {outdir}")
        os.system(f"rm -rf {outdir}")
    if bgzip and not outfile.endswith(ext):
        cmd = f"bgzip {outfile} && tabix -S 1 -s 1 -b 2 -e 3 -f {outfile}.gz"
        print(f"Run bgzip, CMD: {cmd}")
        os.system(cmd)

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
        merge_cz(indir=indir, cz_paths=cz_paths, bgzip=False,
                 outfile=outfile, n_jobs=n_jobs, Path_to_chrom=Path_to_chrom)
# ==========================================================
def extractCG(input=None, outfile=None, ssi=None, chunksize=5000,
              merge_cg=False):
    """
    Extract CG context from .cz file

    Parameters
    ----------
    cz_path :path
    outfile :path
    ssi : path
        ssi should be ssi to mm10_with_chrL.allc.cz.CGN.ssi, not forward
        strand ssi, but after merge (if merge_cg is True), forward ssi
        mm10_with_chrL.allc.cz.+CGN.ssi should be used to generate
         reference, one can
        run: czip extract -m mm10_with_chrL.allc.cz
        -o mm10_with_chrL.allCG.forward.cz
        -b mm10_with_chrL.allc.cz.+CGN.ssi and use
        mm10_with_chrL.allCG.forward.cz as new reference.
    chunksize :int
    merge_cg: bool
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
        if merge_cg:
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


def aggregate(Input=None, Outfile=None, ssi=None, intersect=None, exclude=None,
              chunksize=5000, formats=['H', 'H']):
    """
    Aggregate a given genomic region on a .cz file, for example:
        /usr/bin/time -f "%e\t%M\t%P" czip aggregate -I test.cz -O test_gene.cz \
        -s mm10_with_chrL.allc.genes_flank2k.ssi
    Parameters
    ----------
    Input :
    Outfile :
    ssi :
    intersect :
    exclude :
    chunksize :
    formats :

    Returns
    -------

    """
    cz_path = os.path.abspath(os.path.expanduser(Input))
    ssi_path = os.path.abspath(os.path.expanduser(ssi))
    ssi_reader = Reader(ssi_path)
    reader = Reader(cz_path)
    writer = Writer(Outfile, Formats=formats,
                    Columns=reader.header['Columns'],
                    Dimensions=reader.header['Dimensions'],
                    message=os.path.basename(ssi_path))
    dtfuncs = get_dtfuncs(writer.Formats)
    for dim in reader.dim2chunk_start.keys():
        if dim not in ssi_reader.dim2chunk_start.keys():
            continue
        print(dim)
        IDs = ssi_reader.get_ids_from_ssi(dim)
        # names=[str(record[0], 'utf-8').rstrip('\x00') for record in ssi_reader.__fetch__(dim, s=2, e=3)]
        # if len(IDs.shape) == 1:
        #     records = reader._getRecordsByIds(dim, IDs)
        assert len(IDs.shape) == 2
        records = reader._getRecordsByIdRegions(dim=dim, IDs=IDs)
        data, count = b'', 0
        for record in records:  # unpacked bytes, many values, zip with names
            # record is an array, nrows, two columns (mc and cov)
            sum_v = np.array([0, 0])
            for r in record:  # for every C in a gene region
                sum_v += np.array(struct.unpack(f"<{reader.fmts}", r))
            data += struct.pack(f"<{writer.fmts}",
                                *[func(v) for v, func in zip(sum_v, dtfuncs)])
            count += 1
            if count > chunksize:
                writer.write_chunk(data, dim)
                data, count = b'', 0
        if len(data) > 0:
            writer.write_chunk(data, dim)
    writer.close()
    reader.close()
    ssi_reader.close()

def __split_mat(infile, chrom, snames, outdir, n_ref):
    tbi = pysam.TabixFile(infile)
    records = tbi.fetch(reference=chrom)
    N = n_ref + len(snames) * 2
    fout_dict = {}
    for sname in snames:
        fout_dict[sname] = open(os.path.join(outdir, f"{sname}.{chrom}.bed"), 'w')
        fout_dict[sname].write("chrom\tstart\tend\tstrand\tpval\todd_ratio\n")
    for line in records:
        values = line.replace('\n', '').split('\t')
        if len(values) < N:
            print(infile, chrom)
            raise ValueError("Number of fields is wrong.")
        ch, beg, end, strand = values[:4]
        beg = int(beg)
        end = int(end)
        for i, sname in enumerate(snames):
            or_value = values[n_ref + i * 2]
            try:
                OR = float(or_value)
            except:
                OR = 1
            if OR >= 1:  # hyper methylation
                # continue  # only keep hypomethylated
                pval = 1
            else:
                pval = values[n_ref + i * 2 + 1]
            fout_dict[sname].write(f"{chrom}\t{beg}\t{end}\t{strand}\t{pval}\t{or_value}\n")
    for sname in snames:
        fout_dict[sname].close()
    tbi.close()

def combp(input, outdir="cpv", n_jobs=24, dist=300, temp=True, bed=False):
    """
    Run comb-p on a fisher result matrix (generated by `merge_cz -f fisher`),
    /usr/bin/time -f "%e\t%M\t%P" czip combp -i major_type.fisher.txt.gz -n 64
    Run one samples (all chromosomes), 8308.18(2.3h) 65053112(62G)        321%

    Parameters
    ----------
    input : path
        path to result from merge_cz -f fisher.
    outdir : path
    n_jobs : int
    dist: int
        max distance between two site to be included in one DMR.
    db : str
    temp : bool
        whether to keep temp dir
    bed : bool
        whether to keep bed directory

    Returns
    -------

    """
    try:
        # import cpv
        from cpv.pipeline import pipeline as cpv_pipeline
    except:
        print("Please install cpv using: pip install git+https://github.com/DingWB/combined-pvalues")
    # from multiprocessing import set_start_method
    # set_start_method("spawn")

    infile = os.path.abspath(os.path.expanduser(input))
    outdir = os.path.abspath(os.path.expanduser(outdir))
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    columns = pd.read_csv(infile, sep='\t', nrows=1).columns.tolist()
    snames = [col[:-5] for col in columns[4:] if col.endswith('.pval')]
    tbi = pysam.TabixFile(infile)
    chroms = sorted(tbi.contigs)
    tbi.close()

    bed_dir = os.path.join(outdir, 'bed')
    if not os.path.exists(bed_dir):
        os.mkdir(bed_dir)
        pool = multiprocessing.Pool(n_jobs)
        jobs = []
        print("Splitting matrix into different samples and chroms.")
        for chrom in chroms:
            job = pool.apply_async(__split_mat,
                                   (infile, chrom, snames, bed_dir, 5))
            jobs.append(job)
        for job in jobs:
            r = job.get()
        pool.close()
        pool.join()
    else:
        print("bed directory existed, skip split matrix into bed files.")

    tmpdir = os.path.join(outdir, "tmp")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    pool = multiprocessing.Pool(n_jobs)
    jobs = []
    print("Running cpv..")
    # dist = 750
    acf_dist = int(round(dist / 3, -1))
    # command line, c is 1-based, in API, col_num is 0-based
    for chrom in chroms:
        for sname in snames:
            bed_file = os.path.join(bed_dir, f"{sname}.{chrom}.bed")
            prefix = os.path.join(tmpdir, f"{sname}.{chrom}")
            outfile = os.path.join(tmpdir, f"{sname}.{chrom}.regions-p.bed.gz")
            if os.path.exists(outfile):
                continue
            job = pool.apply_async(cpv_pipeline,
                                   (4, None, dist, acf_dist, prefix,
                                    0.05, 0.05, "refGene", [bed_file], True, 1, None, False,
                                    None, True))
            jobs.append(job)
    for job in jobs:
        r = job.get()
    pool.close()
    pool.join()

    print("Merging cpv results..")
    for sname in snames:
        outfile = os.path.join(outdir, f"{sname}.bed")
        for chrom in chroms:
            infile = os.path.join(tmpdir, f"{sname}.{chrom}.regions-p.bed.gz")
            if not os.path.exists(infile):
                continue
            df = pd.read_csv(infile, sep='\t')
            df = df.loc[df.z_sidak_p <= 0.05]
            if df.shape[0] == 0:
                continue
            if not os.path.exists(outfile):
                df.to_csv(outfile, sep='\t', index=False, header=True)
            else:
                df.to_csv(outfile, sep='\t', index=False, header=False, mode='a')
    merged_dmr_path = os.path.join(outdir, 'merged_dmr.txt')
    data = None
    for sname in snames:
        infile = os.path.join(outdir, f"{sname}.bed")
        df = pd.read_csv(infile, sep='\t')
        df.drop(['min_p', 'z_p', 'z_sidak_p'], inplace=True, axis=1)
        df['sname'] = sname
        if data is None:
            data = df.copy()
        else:
            data = pd.concat([data, df], ignore_index=True)
    # cols = data.columns.tolist()
    # data=data.groupby(cols[:3]).agg(lambda x:x.tolist())
    # data=data.applymap(lambda x:x[0] if len(x)==1 else ','.join(list(map(str,x))))
    data.to_csv(merged_dmr_path, sep='\t', index=False)
    if not bed:
        os.system(f"rm -rf {bed_dir}")
    if not temp:
        os.system(f"rm -rf {tmpdir}")

def annot_dmr(input="merged_dmr.txt", matrix="merged_dmr.cell_class.beta.txt",
              outfile='dmr.annotated.txt', delta_cutoff=None):
    """
    Annotate DMR result from cpv.

    Parameters
    ----------
    dmr : path
        Merged dmr from czip combp.
    matrix :  path
        result of agg_beta using dmr and output of merge_cz (fraction) as input.
    outfile : path
        annotated dmr, containing hypomethylated sname, delta
    Returns
    -------

    """
    data = pd.read_csv(os.path.expanduser(matrix), sep='\t', index_col=[0, 1, 2])
    df_rows = data.index.to_frame()
    # assert data.shape[0] == df_dmr.shape[0]
    # data = data.loc[df_dmr.index.tolist()]
    cols = data.columns.tolist()
    a = data.values
    df_rows['Hypo'] = [cols[i] for i in np.argmin(a, axis=1)]
    df_rows['Hyper'] = [cols[i] for i in np.argmax(a, axis=1)]
    df_rows['Max'] = np.max(a, axis=1)
    df_rows['Min'] = np.min(a, axis=1)
    df_rows['delta_beta'] = df_rows.Max - df_rows.Min
    df_dmr = pd.read_csv(os.path.expanduser(input), sep='\t')
    cols = df_dmr.columns.tolist()
    n_cpg = df_dmr.iloc[:, :4].drop_duplicates().set_index(cols[:3])[cols[3]].to_dict()
    dmr_sample_dict = df_dmr.loc[:, cols[:3] + ['sname']].drop_duplicates().groupby(
        cols[:3]).sname.agg(lambda x: x.tolist())
    df_rows['n_dms'] = df_rows.index.to_series().map(n_cpg)
    df_rows['sname'] = df_rows.index.to_series().map(dmr_sample_dict)
    # df_rows['sname'] = df_rows.apply(lambda x: x.Hypo if x.Hypo in x.sname else np.nan, axis=1)
    # only keep hypomethylated DMR
    # df_rows = df_rows.loc[~ df_rows.sname.isna()]
    df_rows['sname'] = df_rows['sname'].apply(lambda x: ','.join(x))
    if not delta_cutoff is None:
        df_rows = df_rows.loc[df_rows.delta_beta >= delta_cutoff]
    # df_rows.drop(['Hyper'], axis=1, inplace=True)
    df_rows.to_csv(os.path.expanduser(outfile),
                   sep='\t', index=False)

if __name__ == "__main__":
    import fire

    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire()
