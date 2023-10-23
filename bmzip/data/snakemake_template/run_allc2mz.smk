import os.path
import pandas as pd

if 'gcp' in config and config["gcp"]:
    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider

    GS = GSRemoteProvider()
    os.environ[
        'GOOGLE_APPLICATION_CREDENTIALS'] = os.path.expanduser('~/.config/gcloud/application_default_credentials.json')
    print(workflow.default_remote_prefix)

allc_path = config["allc_path"]
allc_files = pd.read_csv(os.path.expanduser(allc_path),sep='\t',
    header=None).iloc[:, 0].tolist()
#allc_path is a file with only one column, no header, if allc files are stored on
# cloud, the final path will be "gs://"+workflow.default_remote_prefix+
# path in allc_path file.
if 'outdir' in config:
    outdir = config["outdir"]  #
else:
    outdir = 'mz'

if not os.path.exists(outdir):
    os.mkdir(outdir)
if 'suffix' in config:
    suffix = config["suffix"]
else:
    # suffix = ".allc.tsv.gz"
    suffix = '.' + '.'.join(os.path.basename(allc_files[0]).split('.')[1:])

if 'indir' in config:
    indir = config["indir"]  #real path = indir + allc_file
else:  #full path
    indir = os.path.dirname(allc_files[0])

snames = [os.path.basename(file.replace(suffix,'')) for file in allc_files]
if 'reference' in config:
    reference = config['reference']
    if 'ref_prefix' in config:
        ref_prefix = config['ref_prefix']
        os.system(f"gsutil cp -n {ref_prefix}/{reference} ./")
else:
    reference = None

print(allc_path,indir,outdir,suffix,reference)

rule target_all:
    input:  #[sname+'.mz' for sname in snames]
        expand(outdir + "/{sample}.mz",sample=snames)

rule run_allc2mz:
    input:
        allc_file=lambda wildcards: os.path.join(indir,wildcards.sname + suffix),
        allc_file_idx=lambda wildcards: os.path.join(indir,wildcards.sname + suffix + '.tbi')

    output:
        "{outdir}/{sname}.mz"
    conda:
        "yap"

    params:
        reference='' if reference is None else f"--reference {reference}",

    shell:
        """
        bmzip allc2mz {input.allc_file} {output} {params.reference}
        """
