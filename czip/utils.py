import os,sys
import struct
from .bmz import Writer,get_dtfuncs
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
if __name__ == "__main__":
    pass
