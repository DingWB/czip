# Installation

# Usage

```shell
/usr/bin/time -f "%e\t%M\t%P" bmzip Writer  -O test.bmzip -F H,H -n mc,cov -t chrom  pack -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 4,5 -c 5000 -t 0

/usr/bin/time -f "%e\t%M\t%P" bmzip Writer -O test_bed.bmzip -F Q,H,H -n pos,mc,cov -t chrom pack -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -c 5000

bmzip summary_blocks
bmzip Reader -f test_bed.bmzip view -t 0
```