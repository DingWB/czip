import os,sys

def extractChr(record, outdir):
	cdef int i,N
	# cdef char* chrom, base, context, strand
	chrom = record.id
	outfile = os.path.join(outdir, chrom + ".bed")
	print(chrom)
	N=record.seq.__len__() - 1
	R=[]
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
		R.append([i + 1,context,strand])
	# position is 0-based (start) 1-based (end position, i+1)
	return chrom,R

def extractC():
