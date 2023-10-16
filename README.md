# Installation
```shell
pip install git+http://github.com/DingWB/bmzip

python setup.py install
```
# Usage

## generate ALLC coordinates
```shell
bmzip AllC -G ~/Ref/mm10/mm10_ucsc_with_chrL.fa -O mm10_with_chrL.mz -j 8 writePattern
```

## Write .mz
```shell
/usr/bin/time -f "%e\t%M\t%P" bmzip Writer  -O test.bmzip -F H,H -N mc,cov -T chrom -v 1 pack -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 4,5 -t 0 -c 5000

/usr/bin/time -f "%e\t%M\t%P" bmzip Writer -O test_bed.bmzip -F Q,H,H -N pos,mc,cov -T chrom pack -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -t 0

bmzip Reader -I test_bed.bmzip summary_blocks
bmzip Reader -I test_bed.bmzip view -t 0
```

### cat multiple .mz
```shell
bmzip Writer -O mm10_with_chrL.mz -F Q,3s,c -N pos,context,strand -T chrom catmz -I "mm10_with_chrL.mz.tmp/*.mz"

bmzip Writer -O mm10_with_chrL.mz -F Q,3s,c -N pos,context,strand -T chrom catmz -I "mm10_with_chrL.mz.tmp/*.mz" -o ~/Ref/mm10/mm10_ucsc_with_chrL.chrom.sizes
```

## View
### print_header
```shell
bmzip Reader -I mm10_with_chrL.mz print_header
# {'magic': b'BMZIP', 'version': 1.0, 'total_size': 3160879498, 'Formats': ['Q', '3s', 'c'], 'names': ['pos', 'context', 'strand'], 'tags': ['chrom'], 'header_size': 51}
```

### view
```shell
bmzip Reader -I test_bed.bmzip view --help
```
```text
INFO: Showing help with the command 'bmzip Reader -I test_bed.bmzip view -- --help'.

NAME
    bmzip Reader -I test_bed.bmzip view - View .bmz file.

SYNOPSIS
    bmzip Reader -I test_bed.bmzip view <flags>

DESCRIPTION
    View .bmz file.

FLAGS
    -s, --show_tag=SHOW_TAG
        Type: Optional[]
        Default: None
        index of tags given to writer.write_chunk, separated by comma, default is None, tags[show_tag] will be shown in each row (such as sampleID and chrom)
    -h, --header=HEADER
        Default: True
        whether to print the header.
    -T, --Tag=TAG
        Type: Optional[]
        Default: None
        None (default): use the default order in .mz file bool (True): sort the tags and used as Tag
```
```shell
bmzip Reader -I mm10_with_chrL.mz view |head
#pos     context strand
#3000003 CTG     +
#3000005 CAG     -
#3000009 CTA     +
#3000016 CAA     -
#3000018 CAC     -
#3000019 CCA     -
#3000023 CTT     +
#3000027 CAA     -
#3000029 CTC     -

# add tags[0] (the first tag) to the output using parameter `--show_tag`
bmzip Reader -I mm10_with_chrL.mz view -s 0 |head
#chrom   pos     context strand
#chr1    3000003 CTG     +
#chr1    3000005 CAG     -
#chr1    3000009 CTA     +
#chr1    3000016 CAA     -
#chr1    3000018 CAC     -
#chr1    3000019 CCA     -
#chr1    3000023 CTT     +
#chr1    3000027 CAA     -
#chr1    3000029 CTC     -

# don't print header
bmzip Reader -I mm10_with_chrL.mz view -s 0 -h False |head 
chr1    3000003 CTG     +
chr1    3000005 CAG     -
chr1    3000009 CTA     +
chr1    3000016 CAA     -

# view .mz data using a different chromosomes order
bmzip Reader -I test_bed.bmzip view -s 0 -h False  |cut -f 1 |uniq
#chr1
#chr10
#chr11
#chr12
#chr13
#chr14
#chr15
#chr16
#chr17
#chr18
#chr19
#chr1_GL456210_random
#chr1_GL456211_random
#chr1_GL456212_random
#chr1_GL456213_random
#chr1_GL456221_random
#chr2
#chr3
#chr4
#chr4_GL456216_random
#chr4_JH584292_random
#chr4_GL456350_random
#chr4_JH584293_random
#chr4_JH584294_random
#chr5
#chr5_JH584296_random
#chr5_JH584297_random
#chr5_JH584298_random
#chr5_GL456354_random
#chr5_JH584299_random
#chr6
#chr7
#chr7_GL456219_random
#chr8
#chr9
#chrM
#chrX
#chrX_GL456233_random
#chrY
#chrY_JH584300_random
#chrY_JH584301_random
#chrY_JH584302_random
#chrY_JH584303_random
#chrUn_GL456239
#chrUn_GL456367
#chrUn_GL456378
#chrUn_GL456381
#chrUn_GL456382
#chrUn_GL456383
#chrUn_GL456385
#chrUn_GL456390
#chrUn_GL456392
#chrUn_GL456393
#chrUn_GL456394
#chrUn_GL456359
#chrUn_GL456360
#chrUn_GL456396
#chrUn_GL456372
#chrUn_GL456387
#chrUn_GL456389
#chrUn_GL456370
#chrUn_GL456379
#chrUn_GL456366
#chrUn_GL456368
#chrUn_JH584304
#chrL

# provide parameter `order` will get a different order
bmzip Reader -I test_bed.bmzip view -s 0 -h False -T ~/Ref/mm10/mm10_ucsc.main.chrom.sizes |cut -f 1 |uniq
#chr1
#chr10
#chr11
#chr12
#chr13
#chr14
#chr15
#chr16
#chr17
#chr18
#chr19
#chr2
#chr3
#chr4
#chr5
#chr6
#chr7
#chr8
#chr9
#chrM
#chrX
#chrY

# only show selected chromosomes using parameter `order`
# this can also be done using query
bmzip Reader -I test_bed.bmzip view -s 0 -h False -T "[['chr1'],['chr2']]" |cut -f 1 |uniq
# chr1
#chr2

# only view 1 chrom
bmzip Reader -I test_bed.bmzip view -s 0 -h False -T chr1 |head
```

