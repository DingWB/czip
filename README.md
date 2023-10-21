# Block Methylation ZIP

## Installation
```shell
pip install git+http://github.com/DingWB/bmzip

python setup.py install

conda install -c "intel/label/test" libfabric ?
```

## Implementation
![docs/images/img.png](docs/images/img.png)

## Usage

### generate ALLC coordinates (.mz file was created by ALLC class)
```shell
bmzip AllC -G ~/Ref/mm10/mm10_ucsc_with_chrL.fa -O mm10_with_chrL.allc.mz -j 20 run
# took 15 minutes using 20 cpus
```

### Create .mz using `bmzip Writer`
```shell
/usr/bin/time -f "%e\t%M\t%P" bmzip Writer  -O test.mz -F H,H -C mc,cov -D chrom -v 1 tomz -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 4,5 -d 0 -c 5000

# pack from stdin
zcat FC_E17a_3C_1-1-I3-F13.allc.tsv.gz |cut -f 1,5,6 | bmzip Writer -O test.mz -F H,H -C mc,cov -D chrom -v 1 tomz -I stdin -u 1,2 -d 0

/usr/bin/time -f "%e\t%M\t%P" bmzip Writer -O test_bed.mz -F Q,H,H -C pos,mc,cov -D chrom tomz -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -d 0

bmzip Reader -I test_bed.mz summary_blocks
bmzip Reader -I test_bed.mz view -s 0
```

#### cat multiple .mz files into one .mz file
```shell
bmzip Writer -O mm10_with_chrL.allc.mz -F Q,c,3s -C pos,strand,context -D chrom catmz -I "mm10_with_chrL.allc.mz.tmp/*.mz"

bmzip Writer -O mm10_with_chrL.allc.mz -F Q,c,3s -C pos,strand,context -D chrom catmz -I "mm10_with_chrL.allc.mz.tmp/*.mz" --dim_order ~/Ref/mm10/mm10_ucsc_with_chrL.chrom.sizes
```

### Using reference coordinates when create .mz file (without coordinates)
For single cell DNA methylation datasets, we can create .mz files only contain mv and cov, no coordinates, cause all cells share the same set of coordinates, which can be created using  bmzip AllC.
```shell
bmzip Writer  -O test.mz -F H,H -C mc,cov -D chrom -v 1 tomz -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u '[4,5]' -d '[0]' -r mm10_with_chrL.allc.mz -pr "['pos']" -p '[1]' -j 8
```

### convert allc to .mz
```shell
# with reference
bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz test.mz -r mm10_with_chrL.allc.mz -v 1

# without reference
bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz test.mz -v 1

# with coordinates
bmzip Writer -O test_bed.mz -F Q,H,H -C pos,mc,cov -D chrom tomz -I FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -d 0

# allc2mz in parallel
time bmzip allc2mz -r mm10_with_chrL.allc.mz -n 48 -P ~/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt allc_path.txt test
```

### test difference

```shell
# Test difference between original allc.tsv.gz against .mz with position
bmzip Writer  -O test_bed1.mz -F Q,H,H -C pos,mc,cov -D chrom -v 1 tomz -I FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -d 0
bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz test_bed2.mz -v 1
bmzip test_diff --file1 FC_E17a_3C_1-1-I3-F13.allc.tsv.gz --usecols1 "[0,1,4,5]" --header1 0 --file2 test_bed1.tsv --usecols2 "[0,1,2,3]" --header2 1
bmzip test_diff --file1 FC_E17a_3C_1-1-I3-F13.allc.tsv.gz --usecols1 "[0,1,4,5]" --header1 0 --file2 test_bed2.tsv --usecols2 "[0,1,2,3]" --header2 1

# Test difference between original allc.tsv.gz against .mz without position (using reference)
bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz test_bed3.mz -v 1 -r mm10_with_chrL.allc.mz

# pack from stdin
zcat FC_E17a_3C_1-1-I3-F13.allc.tsv.gz |cut -f 1,5,6 | bmzip Writer -O test.mz -F H,H -C mc,cov -D chrom -v 1 tomz -I stdin -u 1,2 -d 0

bmzip Writer -O test_bed.mz -F Q,H,H -C pos,mc,cov -D chrom tomz -I /anvil/scratch/x-wding2/Projects/mouse-pfc/test_ballc/test_bmzip/FC_E17a_3C_1-1-I3-F13.allc.tsv.gz -u 1,4,5 -d 0
```

### View

#### print_header

```shell
bmzip Reader -I mm10_with_chrL.allc.mz print_header
# {'magic': b'BMZIP', 'version': 1.0, 'total_size': 3160879498, 'Formats': ['Q', '3s', 'c'], 'names': ['pos', 'context', 'strand'], 'tags': ['chrom'], 'header_size': 51}
```

#### view

```shell
bmzip Reader -I test_bed.mz view --help
```
```text
INFO: Showing help with the command 'bmzip Reader -I test_bed.mz view -- --help'.

NAME
    bmzip Reader -I test_bed.mz view - View .bmz file.

SYNOPSIS
    bmzip Reader -I test_bed.mz view <flags>

DESCRIPTION
    View .bmz file.

FLAGS
    -s, --show_dim=SHOW_DIM
        Type: Optional[]
        Default: None
        index of dims given to writer.write_chunk, separated by comma, default is None, dims[show_dim] will be shown in each row (such as sampleID and chrom)
    -h, --header=HEADER
        Default: True
        whether to print the header.
    -d, --dim=DIM
        Type: Optional[]
        Default: None
        None (default): use the default order in .mz file;

        bool (True): sort the dims and used as dim
```
```shell
bmzip Reader -I mm10_with_chrL.allc.mz view |head
#pos     strand  context
#3000003 +       CTG
#3000005 -       CAG
#3000009 +       CTA
#3000016 -       CAA
#3000018 -       CAC

# add tags[0] (the first tag) to the output using parameter `--show_tag`
bmzip Reader -I mm10_with_chrL.allc.mz view -s 0 |head
#chrom   pos     context strand
#chr1    3000003 CTG     +
#chr1    3000005 CAG     -
#chr1    3000009 CTA     +
#chr1    3000016 CAA     -

# don't print header
bmzip Reader -I mm10_with_chrL.allc.mz view -s 0 -h False |head 
chr1    3000003 CTG     +
chr1    3000005 CAG     -
chr1    3000009 CTA     +

# view .mz data using a different chromosomes order
# default tags order
bmzip Reader -I test_bed.mz view -s 0 -h False  |cut -f 1 |uniq
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
#chr1_GL456221_random
#chr2
#chr3
#chr4
#chr4_GL456216_random
#chr4_JH584294_random
#chr5
#chr5_JH584296_random
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
#chrY_JH584303_random
#chrUn_GL456239
#chrUn_GL456385
#chrL

# provide parameter `dim` will get a different order
# use the dim (chrom) order from chrom size file (first columns).
bmzip Reader -I test_bed.mz view -s 0 -h False -d ~/Ref/mm10/mm10_ucsc.main.chrom.sizes |cut -f 1 |uniq
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

# only show selected chromosomes using parameter `dim`
# this can also be done using query
bmzip Reader -I test_bed.mz view -s 0 -h False -d chr1,chr2 |cut -f 1 |uniq
# chr1
#chr2

# only view 1 chrom
bmzip Reader -I test_bed.mz view -s 0 -h False -d chr1 |head
```

### Query
```shell
/usr/bin/time -f "%e\t%M\t%P" bmzip Reader -I test_bed.mz query -d chr8 -s 129300305 -e 129300362
#chrom   pos     mc      cov
#chr8    129300305       0       1
#chr8    129300317       0       1
#chr8    129300330       0       1
#chr8    129300333       0       1
#chr8    129300339       0       1
#chr8    129300347       0       1
#chr8    129300352       0       1
#chr8    129300355       0       1
#chr8    129300362       0       1

#0.80    239428  118%
# 129300305-129300362 is at the last position of chr8.


# compare with the original allc file
zcat FC_E17a_3C_1-1-I3-F13.allc.tsv.gz |awk '$5 >=100' |head
#chr12   3109883 +       CGG     224     255     1
#chr12   3109884 -       CGT     590     690     1
#chr12   3109911 +       CGT     410     471     1
#chr12   3109912 -       CGT     1364    1548    1
#chr12   3109923 +       CGA     451     505     1
#chr12   3109924 -       CGC     1413    1609    1
#chr12   3109971 +       CAT     398     522     1
#chr12   3110026 +       CGA     401     457     1
#chr12   3110027 -       CGT     1218    1356    1
#chr12   3110041 +       CAA     263     398     1

/usr/bin/time -f "%e\t%M\t%P" bmzip Reader -I test_bed.mz query -D chr12 -s 3109883 -e 3110041
#chrom   pos     mc      cov
#chr12   3109883 224     255
#chr12   3109884 590     690
#chr12   3109885 20      815
#chr12   3109891 2       939
#chr12   3109893 2       959
#chr12   3109899 9       404
#chr12   3109901 8       420
#chr12   3109903 8       430
#chr12   3109908 15      1437
#chr12   3109909 11      1453
#chr12   3109911 410     471
#chr12   3109912 1364    1548
#chr12   3109914 10      1579
#chr12   3109921 19      1610
#chr12   3109922 4       1513
#chr12   3109923 451     505
#chr12   3109924 1413    1609
#chr12   3109926 8       1631
#chr12   3109927 9       1623
#chr12   3109932 3       537
#chr12   3109934 14      1626
#chr12   3109940 9       1628
#chr12   3109941 8       1585
#chr12   3109943 19      1637
#chr12   3109944 9       1649
#chr12   3109953 13      1647
#chr12   3109958 14      1654
#chr12   3109960 2       527
#chr12   3109961 7       530
#chr12   3109963 3       531
#chr12   3109965 21      1603
#chr12   3109968 11      1645
#chr12   3109969 14      1643
#chr12   3109971 398     522
#chr12   3109974 21      1638
#chr12   3109975 14      1639
#chr12   3109981 28      1635
#chr12   3109982 18      1625
#chr12   3109983 6       529
#chr12   3109986 11      1635
#chr12   3109991 2       517
#chr12   3109993 11      1639
#chr12   3109999 9       522
#chr12   3110002 26      1624
#chr12   3110003 15      1617
#chr12   3110009 17      1599
#chr12   3110011 8       1571
#chr12   3110015 3       477
#chr12   3110018 4       475
#chr12   3110019 4       484
#chr12   3110021 5       470
#chr12   3110024 18      1427
#chr12   3110026 401     457
#chr12   3110027 1218    1356
#chr12   3110029 3       443
#chr12   3110032 13      1224
#chr12   3110039 5       990
#chr12   3110041 263     398
#0.98    238200  90%, only took 0.98 s


bmzip Reader -I 2.mz query -D "{'chrom':'chr12'}" -s 3110029 -e 3110041 -r mm10_with_chrL.allc.mz
```

#### Query with reference

```shell
bmzip1 Reader -I test_bed3.mz query -d chr12 -s 3109911 -e 3109913 -r mm10_with_chrL.allc.mz |head

bmzip1 Reader -I mm10_with_chrL.allc.mz query -d chr12 -s 3109911 -e 3109913  |head
```

### SubSet Index (.ssi)

```shell
bmzip generate_context_ssi -I mm10_with_chrL.allc.mz -p CGN
```