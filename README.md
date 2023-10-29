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
bmzip AllC -G ~/Ref/mm10/mm10_ucsc_with_chrL.fa -O mm10_with_chrL.allc.mz -n 20 -k True run
# took 15 minutes using 20 cpus

bmzip generate_ssi ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz
# output is ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz.+CGN.bmi

bmzip extract -m ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz -o ~/Ref/mm10/annotations/mm10_with_chrL.allCG.forward.mz -b ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz.+CGN.bmi
#output is: mm10_with_chrL.allCG.forward.mz
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
time bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz FC_E17a_3C_1-1-I3-F13.mz -r mm10_with_chrL.allc.mz
time bmzip allc2mz FC_P13a_3C_4-4-C7-E19.allc.tsv.gz FC_P13a_3C_4-4-C7-E19.mz -r mm10_with_chrL.allc.mz

# without reference
bmzip allc2mz FC_E17a_3C_1-1-I3-F13.allc.tsv.gz test.mz -v 1

# allc2mz in parallel
time bmzip allc2mz -r mm10_with_chrL.allc.mz -n 48 -P ~/Ref/mm10/mm10_ucsc_with_chrL.main.chrom.sizes.txt allc_path.txt test

# make a small example allc file in which there is only 1 record in each chrom
zcat FC_E17a_3C_1-1-I3-F13.allc.tsv.gz |awk 'BEGIN{chrom=""};{if ($1!=chrom) {print($0); chrom=$1}}' | bgzip > 1.allc.tsv.gz
tabix -b 2 -e 2 1.allc.tsv.gz
time bmzip allc2mz 1.allc.tsv.gz 1.mz -r mm10_with_chrL.allc.mz

# file sizes
24514 allc.tsg.gz + .tbi
2.83TiB

.mz files:
769.21GiB
```

### Extract CG from .mz and merge strand

```shell
bmzip extractCG -I mz/FC_P13a_3C_2-1-E5-D13.mz -o FC_P13a_3C_2-1-E5-D13.CGN.mz -b ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz.CGN.bmi

# create reference for forward CGN
bmzip extract -m ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz -o ~/Ref/mm10/annotations/mm10_with_chrL.allCG.forward.mz -b ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz.+CGN.bmi -o mm10_with_chrL.allCG.forward.mz

# view CG .mz
bmzip Reader -I FC_P13a_3C_2-1-E5-D13.CGN.mz view -s 0 -r ~/Ref/mm10/annotations/mm10_with_chrL.allCG.forward.mz
```

### Run allc2mz using snakemake

```shell
# local
snakemake -s /home/x-wding2/Projects/Github/bmzip/data/snakemake_template/run_allc2mz.smk --config allc_path="allc_path.txt" reference="mm10_with_chrL.allc.mz" -j 2 -np

# gcp
head 1000_allc_path.tsv #gs://mouse_pfc/allc/FC_P0a_3C_10-6-G10-D11.allc.tsv.gz
#FC_P0a_3C_10-6-G10-D11.allc.tsv.gz
#FC_P0a_3C_8-2-I2-P15.allc.tsv.gz

snakemake --printshellcmds --immediate-submit -s run_allc2mz.smk --config indir="allc" outdir="mz" allc_path="1000_allc_path.tsv" reference="mm10_with_chrL.allc.mz" ref_prefix="gs://wubin_ref/mm10/annotations" Path_to_chrom="mm10_ucsc_with_chrL.main.chrom.sizes.txt" Path_to_chrom_prefix="gs://wubin_ref/mm10" gcp=True --default-remote-prefix mouse_pfc --default-remote-provider GS --google-lifesciences-region us-west1 --scheduler greedy -j 96 -np
# --keep-remote

bmzip copy_smk -o allc2mz.smk
snakemake --printshellcmds --immediate-submit --notemp -s allc2mz.smk --config indir="gs://mouse_pfc/test_allc" outdir="test_mz" reference="mm10_with_chrL.allc.mz" ref_prefix="gs://wubin_ref/mm10/annotations" chrom="mm10_ucsc_with_chrL.main.chrom.sizes.txt" chrom_prefix="gs://wubin_ref/mm10" gcp=True --default-remote-prefix mouse_pfc --default-remote-provider GS --google-lifesciences-region us-west1 --scheduler greedy -j 96 -np

# indir is either a GCP cloud dirname(gs://mouse_pfc/allc, no need for allc_path) or a dir prefix (such as allc, used with allc_path)
bmzip prepare_sky --indir allc --outdir pfc_mz --allc_path allc.path --reference mm10_with_chrL.allc.mz --ref_prefix gs://wubin_ref/mm10/annotations --gcp True --bucket mouse_pfc --chrom mm10_ucsc_with_chrL.main.chrom.sizes.txt --chrom_prefix gs://wubin_ref/mm10 --cpu 96 > 1.yaml
sky spot launch -n allc2mz 1.yaml
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

# examine total number of rows
bmzip Reader -I mm10_with_chrL.allc.mz summary_chunks |awk 'BEGIN{sum=0};{if(NR > 1){sum+=$6}};END{print(sum)}'
# should be 1105362117 for mm10
```

### Query
```shell
bmzip Reader -I test.mz query -D chr8 -s 129300305 -e 129300362
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

bmzip Reader -I test_no_ref.mz query -D chr12 -s 3109883 -e 3110041
#chrom   pos     mc      cov
#chr12   3109883 224     255
#chr12   3109884 590     690
#chr12   3109885 20      815
#chr12   3109991 2       517
#chr12   3109993 11      1639
#chr12   3109999 9       522
#chr12   3110002 26      1624
#chr12   3110003 15      1617
#chr12   3110009 17      1599

#0.98    238200  90%, only took 0.98 s


bmzip Reader -I FC_E17a_3C_1-1-I3-F13.mz query -D "{'chrom':'chr12'}" -s 3110029 -e 3110041 -r mm10_with_chrL.allc.mz
tabix FC_E17a_3C_1-1-I3-F13.allc.tsv.gz chr12:3110029-3110041

bmzip Reader -I test.mz query -D chr11 -s 3104261 -e 3104273 -r mm10_with_chrL.allc.mz
tabix FC_E17a_3C_1-1-I3-F13.allc.tsv.gz chr11:3104261-3104273

tabix FC_P13a_3C_4-4-C7-E19.allc.tsv.gz chr10:74629843-74629846
bmzip Reader -I test_small_file.mz query -D "{'chrom':'chr10'}" -s 74629843 -e 74629846 -r mm10_with_chrL.allc.mz
```

#### Query with reference

```shell
bmzip Reader -I test.mz query -D chr12 -s 3109911 -e 3109913 -r mm10_with_chrL.allc.mz |head

bmzip Reader -I mm10_with_chrL.allc.mz query -D chr12 -s 3109911 -e 3109913  |head
```

### SubSet

#### Create subset index (.ssi)

```shell
bmzip generate_context_ssi -I ~/Ref/mm10/annotations/mm10_with_chrL.allc.mz -p CGN
```

#### view subset

```shell
bmzip Reader -I FC_E17a_3C_1-1-I3-F13.mz subset -d chr1 -b mm10_with_chrL.allc.mz.CGN.bmi -r mm10_with_chrL.allc.mz |head

zcat FC_E17a_3C_1-1-I3-F13.allc.tsv.gz |awk '$4 ~ "^CG" && $5 > 3' |head

```

### Aggregration

#### catmz

```shell
bmzip1 Writer -O two_samples.mz -F H,H -C mc,cov -D chrom,filename catmz -I "raw/*.mz" -a
bmzip1 Reader -I two_samples.mz summary_chunks
```