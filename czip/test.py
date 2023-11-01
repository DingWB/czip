from czip import bed2cz
import cProfile


def test_bed2cz(input=None, outfile="test.cz", ref=None, chrom=None):
    if input is None:
        input = "/home/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz"
    if ref is None:
        ref = "/home/x-wding2/Ref/mm10/annotations/mm10_with_chrL.allc.cz"
    if chrom is None:
        chrom = "/home/x-wding2/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt"
    cProfile.run(f"bed2cz('{input}', '{outfile}', reference='{ref}', missing_value=[0, 0],\
           Formats=['H', 'H'], Columns=['mc', 'cov'], Dimensions=['chrom'],\
           usecols=[4, 5], pr=0, pa=1, sep='\t', Path_to_chrom='{chrom}')")


if __name__ == "__main__":
    import fire

    fire.core.Display = lambda lines, out: print(*lines, file=out)
    fire.Fire()
