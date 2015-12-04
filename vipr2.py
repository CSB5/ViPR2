#!/usr/bin/env python
"""FIXME:add-doc

"""


import os
import sys
import argparse
import logging
import copy
import json
import shutil
import subprocess
import getpass

try:
    # python 3
    from itertools import zip_longest
except ImportError:
    # python 2
    from itertools import izip_longest as zip_longest


__author__ = "Andreas Wilm"
__version__ = "2.0.0a"
__email__ = "wilma@gis.a-star.edu.sg"
__license__ = "The MIT License (MIT)"


# global logger
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
    format='%(levelname)s [%(asctime)s]: %(message)s')


# The wrapper script created here calling snakemake
SNAKEMAKE_CLUSTER_WRAPPER = "snake.sh"
# root name of snakemake's make file
SNAKEMAKE_FILE = "snake.make"
# path to snakemake's makefile template
SNAKEMAKE_TEMPLATE = os.path.join(
    os.path.dirname(sys.argv[0]), SNAKEMAKE_FILE)
# root name of config file loaded in SNAKEMAKE_FILE


CONFIG_FILE = "conf.json"
# template variables to be written to CONFIG_FILE
# FIXME might be better to have a template as well and add user variables only
CONF = dict()
# databases and files
CONF['FAMAS'] = "/mnt/software/stow/famas-0.0.10/bin/famas"
CONF['BWA'] = '/mnt/software/stow/bwa-0.7.12/bin/bwa'
CONF['DEBUG'] = False
# CONF['REFFA'] written later
# CONF['SAMPLENAME'] written later
# CONF['SAMPLES'] written later

def main():
    """The main function
    """

    for f in CONF.keys():
        if f in ['DEBUG']:
            continue
        if not os.path.exists(CONF[f]):
            LOG.fatal("Missing file: {}".format(CONF[f]))
            sys.exit(1)
    assert os.path.exists(SNAKEMAKE_TEMPLATE)


    parser = argparse.ArgumentParser(description='VIPR: version 2')
    parser.add_argument('-1', "--fq1", required=True, nargs="+",
                        help="Paired-end FastQ file #1 (gzip supported). Multiple (split) input files allowed")
    parser.add_argument('-2', "--fq2", required=True, nargs="+",
                        help="Paired-end FastQ file #2 (gzip supported). Multiple (split) input files allowed")
    parser.add_argument('-o', "--outdir", required=True,
                        help='Output directory (may not exist, unless using --continue)')
    parser.add_argument('-r', "--reffa", required=True,
                        help='Reference genome')
    parser.add_argument('-n', "--name", required=True,
                        help='Sample name (used as name for assembled genome)')
    default = 8
    parser.add_argument('-c', '--num-cores', type=int, default=default,
                        help='Number of cores to use (default = {})'.format(default))
    parser.add_argument('--verbose', action="store_true",
                        help='Be verbose')
    parser.add_argument('--debug', action="store_true",
                        help='Enable debugging output')
    parser.add_argument('--no-run', action="store_true",
                        help="Prepare output directory and files but don't actually submit job")
    args = parser.parse_args()


    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    if os.path.exists(args.outdir):
        LOG.fatal("Output directory must not exist: {}".format(args.outdir))
        sys.exit(1)

    assert os.path.exists(args.reffa)

    if args.fq1 == args.fq2:
        LOG.fatal("Paired-End FastQ files have identical names")
        sys.exit(1)
    for fq1, fq2 in zip_longest(args.fq1, args.fq2):
        # only i|zip_longest uses None if one is missing
        if fq1 is None or fq2 is None:
            LOG.fatal("Unequal number of FastQ files")
            sys.exit(1)
        # enforce gzipped fastq
        # check they all exist
        for f in [fq1, fq2]:
            if not os.path.exists(f):
                LOG.fatal("FastQ file {} does not exist".format(f))
                sys.exit(1)
            elif not f.endswith(".gz"):
                LOG.fatal("Non-gzipped FastQ files not supported")
                sys.exit(1)

    os.mkdir(args.outdir)

    samples = []
    fqs1 = [os.path.abspath(f) for f in args.fq1]
    fqs2 = [os.path.abspath(f) for f in args.fq2]
    #if not args.no_fastq_sort:
    fqs1 = sorted(fqs1)
    fqs2 = sorted(fqs2)
    for i, (fq1, fq2) in enumerate(zip(fqs1, fqs2)):
        os.symlink(os.path.abspath(fq1), os.path.join(args.outdir, "{}_R1.fastq.gz".format(i+1)))
        os.symlink(os.path.abspath(fq2), os.path.join(args.outdir, "{}_R2.fastq.gz".format(i+1)))
        samples.append("{}_".format(i+1))


    config_file = os.path.join(args.outdir, CONFIG_FILE)
    conf = copy.deepcopy(CONF)
    conf['SAMPLES'] = samples
    conf['REFFA'] = args.reffa
    conf['SAMPLENAME'] = args.samplename.replace(" ", "_")

    with open(config_file, 'w') as fh:
        json.dump(conf, fh, indent=4)

    shutil.copyfile(SNAKEMAKE_TEMPLATE, 
                    os.path.join(args.outdir, SNAKEMAKE_FILE))

    snakemake_cluster_wrapper = os.path.join(args.outdir, SNAKEMAKE_CLUSTER_WRAPPER)
    mail_option = "-m bes -M {}@gis.a-star.edu.sg".format(getpass.getuser())
    with open(snakemake_cluster_wrapper, 'w') as fh:
        fh.write('# snakemake requires python3\n')
        fh.write('source activate py3k;\n')
        fh.write('cd {};\n'.format(os.path.abspath(args.outdir)))
        fh.write('# qsub for snakemake itself\n')
        fh.write('qsub="qsub -pe OpenMP 1 -l mem_free=1G -l h_rt=48:00:00 {} -j y -V -b y -cwd";\n'.format(mail_option))
        fh.write('# -j in cluster mode is the maximum number of spawned jobs\n')
        fh.write('$qsub -N snakemake -o snakemake.qsub.log')
        qsub_per_task = "qsub -pe OpenMP {threads} -l mem_free=8G -l h_rt=24:00:00 -j y -V -b y -cwd"
        fh.write(' \'snakemake -j 8 -c "{}" -s {} --configfile {}\';\n'.format(
               qsub_per_task, SNAKEMAKE_FILE, CONFIG_FILE))
        # FIXME max runtime and mem should be defined per target in SNAKEMAKE_FILE

    cmd = ['bash', snakemake_cluster_wrapper]
    if args.no_run:
        print("Not actually submitting job")
        print("When ready run: {}".format(' '.join(cmd)))
    else:
        print("Running: {}".format(cmd))
        subprocess.check_output(cmd)
                

if __name__ == '__main__':
    main()

