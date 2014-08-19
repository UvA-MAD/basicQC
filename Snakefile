import pandas as pd

SFF_STORE = "/zfs/datastore0/group_root/MAD-RBAB/04_MAD-RBAB-runs/data"
ANALYSIS_HOME = "./"
DESIGN_FILE = "seqdesign.txt"

#  get sample names from design file
df = pd.read_csv(DESIGN_FILE, sep="\t")
colnames = list(df.columns)
colnames[0] = colnames[0].lower()
df.columns = colnames
SAMPLES = list(set(df['sampleid']))

rule get_fq:
    input: DESIGN_FILE
    output: "./raw/{s}.fastq".format(s=s) for s in SAMPLES
    shell: "python get_experiment_fq.py -d {DESIGN_FILE} -s {SFF_STORE} -a {ANALYSIS_HOME} -p"
