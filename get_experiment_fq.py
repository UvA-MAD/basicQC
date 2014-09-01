import subprocess
import argparse
import pandas as pd
import re
import os
import sys


def sample_sff2fq(sample_info):
    """
    For each sampleid generate fatq file(s) from sff.
    """
    # unpack arguments
    sample_data, sff_home, analysis_home = sample_info
    raw_dir = os.path.join(analysis_home, "raw")
    # patttern of sff file names
    sff_pattern = '^IonXpress(RNA)?_0%s_R_[0-9]{4}_([0-9]{2}_){5}user_PRO-[0-9]*-%s.sff$'
    sampleid = list(sample_data['sampleid'])[0]
    sample_sffs = []
    for index, row in sample_data.iterrows():
        rid = row['runid']
        bc = row['barcode'][2:]
        subdir, dirs, files = next(os.walk(os.path.join(sff_home, rid)))

        sff_regexp = re.compile(sff_pattern % (bc, rid))
        for f in files:
            if sff_regexp.match(f):
                sff_file = os.path.join(sff_home, rid, f)
                sample_sffs.append(sff_file)

    # generate fastq files from sff files for this sample
    sample_fastqs = []
    for i, sff_file in enumerate(sample_sffs):
        fastq_file = os.path.join(
            analysis_home,
            '_'.join((str(sampleid), str(i))) + '.fastq'
            )
        sample_fastqs.append(fastq_file)
        subprocess.call(["sff2fastq", sff_file, "-o", fastq_file])
    # concatenate fastq for same sample from multiple runs
    sample_fastq = os.path.join(raw_dir, str(sampleid) + '.fastq')
    cat_command = ['cat'] + sample_fastqs
    with open(sample_fastq, 'w') as fout:
        subprocess.call(cat_command, stdout=fout)

    # remove fastq  files from before concatenation
    [os.remove(f) for f in sample_fastqs]


def main(argv):
    # setup argument parse
    parser = argparse.ArgumentParser(
        description="Generate fatq files from sff files.",
        epilog="Finds sff files for samples described in design file."
               "If sample has been run in muliple runs, fq files are combine to one."
               "Generated fastq files are placed in ./raw relative to 'analysis_home'.")

    parser.add_argument('-d', '--design_file',
                        type=argparse.FileType('r'),
                        help="path to design file")
    parser.add_argument('-s', '--sff_home',
                        type=str,
                        help="path to original sff files")
    parser.add_argument('-a', '--analysis_home',
                        type=str,
                        help="path to target directory with analysis")
    parser.add_argument('-n', '--numthreads',
                        default=1,
                        type=int,
                        help="number of threads, max 10")

    # print help if no command line arguments have been used
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    # parse arguments
    args = parser.parse_args(argv)

    # parse design csv, which is tab separated file
    df = pd.read_csv(args.design_file, sep="\t")

    # first three column names should be case insensitive
    # therefore convert to lower case
    colnames = list(df.columns)
    colnames[:3] = [s.lower() for s in colnames[:3]]
    df.columns = colnames

    # get unique sample names
    # sample name can repeat it design file
    samples = list(set(df['sampleid']))

    # create dir for fastq files
    raw_dir = os.path.join(args.analysis_home, "raw")
    if not os.path.exists(raw_dir):
        os.makedirs(raw_dir)

    numthreads = int(args.numthreads)
    if numthreads > 1:
        # max threads is 10
        n = (10 if numthreads > 10 else numthreads)
        import multiprocessing
        arguments = []
        for sampleid in samples:
            # rows corresponding to particular sample
            # more than one if sample has been run in multiple runs
            sample_data = df[df['sampleid'] == sampleid]
            arguments.append((sample_data,
                              args.sff_home,
                              args.analysis_home))

        pool = multiprocessing.Pool(processes=n)
        pool.map(sample_sff2fq, arguments)
    else:
        for sampleid in samples:
            sample_data = df[df['sampleid'] == sampleid]
            sample_sff2fq((sample_data, args.sff_home, args.analysis_home))

if __name__ == "__main__":
    main(sys.argv[1:])
