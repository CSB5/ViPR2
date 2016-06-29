#!/usr/bin/env python3
"""Successor of the Viral Pipeline Runner (ViPR; see https://github.com/CSB5/vipr)

Assembles your viral amplicons sequences, maps reads and calls low frequency variants
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

LOG_REL_DIR = "logs"


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
logger = logging.getLogger("")
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
# programs
CONF['FAMAS'] = "/mnt/software/stow/famas-0.0.11/bin/famas"
CONF['BWA'] = '/mnt/software/stow/bwa-0.7.12/bin/bwa'
CONF['SAMTOOLS'] = '/mnt/software/stow/samtools-1.3/bin/samtools'
CONF['IVA'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "iva.sh"))
CONF['PRIMER_POS_FROM_SEQ'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "primer_pos_from_seq.sh"))
CONF['COVERAGE_PLOT'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "coverage_plot.py"))
#CONF['MAPPING_SUCCESS'] = os.path.abspath(
#        os.path.join(os.path.dirname(sys.argv[0]), "mapping_success.sh"))
CONF['PRIMER_POS_FROM_SEQ'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "primer_pos_from_seq.sh"))
CONF['MARK_PRIMER'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "mark_primer.py"))
CONF['PRIMER_POS_TO_BED'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "primer_pos_to_bed.py"))
CONF['LOFREQ'] = "/mnt/software/stow/lofreq_star-2.1.2/bin/lofreq"
CONF['BAMLEFTALIGN'] = "/mnt/software/stow/freebayes-1.0.1/bin/bamleftalign"
CONF['SIMPLE_CONTIG_JOINER'] = "/mnt/software/stow/simple-contig-joiner-0.3/bin/simple_contig_joiner.py"

CONF['PRIMER_LEN'] = 25
CONF['VCF2CSV'] = os.path.abspath(
    os.path.join(os.path.dirname(sys.argv[0]), "vcf2csv.py"))
CONF['MUMMERDIR'] = "/mnt/software/unstowable/mummer-3.23/"

# settings
CONF['DEBUG'] = False
# written dynamically
# CONF['REFFA']
# CONF['SAMPLENAME']
# CONF['SAMPLES']
# config[PRIMER_FILE]


def main():
    """The main function
    """

    for f in CONF.keys():
        if f in ['DEBUG', 'PRIMER_LEN']:
            continue
        if not os.path.exists(CONF[f]):
            logger.fatal("Missing file: %s", CONF[f])
            sys.exit(1)
    assert os.path.exists(SNAKEMAKE_TEMPLATE)


    parser = argparse.ArgumentParser(description='VIPR: version 2')
    parser.add_argument('-1', "--fq1", required=True, nargs="+",
                        help="Paired-end FastQ file #1 (gzip only). Multiple (split) input files allowed")
    parser.add_argument('-2', "--fq2", required=True, nargs="+",
                        help="Paired-end FastQ file #2 (gzip only). Multiple (split) input files allowed")
    parser.add_argument('-o', "--outdir", required=True,
                        help='Output directory (may not exist)')
    parser.add_argument('-r', "--reffa", required=True,
                        help='Reference genome')
    parser.add_argument('-p', "--primers", required=True,
                        help='Fasta file containing primers')
    parser.add_argument('-n', "--name", required=True,
                        help='Sample name (used file and sequence naming)')
    parser.add_argument('--verbose', action="store_true",
                        help='Be verbose')
    parser.add_argument('--debug', action="store_true",
                        help='Enable debugging output')
    parser.add_argument('--no-run', action="store_true",
                        help="Prepare output directory and files but don't actually submit job")
    args = parser.parse_args()


    if args.verbose:
        logger.setLevel(logging.INFO)
    if args.debug:
        logger.setLevel(logging.DEBUG)

    if os.path.exists(args.outdir):
        logger.fatal("Output directory must not exist: %s", args.outdir)
        sys.exit(1)

    for f in [args.reffa, args.primers]:
        if not os.path.exists(f):
            logger.fatal("Missing input file %s", f)
            sys.exit(1)

    if args.fq1 == args.fq2:
        logger.fatal("Paired-End FastQ files have identical names")
        sys.exit(1)
    for fq1, fq2 in zip_longest(args.fq1, args.fq2):
        # only i|zip_longest uses None if one is missing
        if fq1 is None or fq2 is None:
            logger.fatal("Unequal number of FastQ files")
            sys.exit(1)
        # enforce gzipped fastq
        # check they all exist
        for f in [fq1, fq2]:
            if not os.path.exists(f):
                logger.fatal("FastQ file %s does not exist", f)
                sys.exit(1)
            elif not f.endswith(".gz"):
                logger.fatal("Non-gzipped FastQ files not supported")
                sys.exit(1)

    # create outdir and logs
    os.makedirs(os.path.join(args.outdir, LOG_REL_DIR))

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
    conf['REFFA'] = os.path.abspath(args.reffa)
    conf['SAMPLENAME'] = args.name.replace(" ", "_")
    conf['PRIMER_FILE'] = os.path.abspath(args.primers)# FIXME copy?

    with open(config_file, 'w') as fh:
        json.dump(conf, fh, indent=4)

    shutil.copyfile(SNAKEMAKE_TEMPLATE,
                    os.path.join(args.outdir, SNAKEMAKE_FILE))

    snakemake_cluster_wrapper = os.path.join(args.outdir, SNAKEMAKE_CLUSTER_WRAPPER)
    mail_option = "-m bes -M {}@gis.a-star.edu.sg".format(getpass.getuser())
    with open(snakemake_cluster_wrapper, 'w') as fh:
        fh.write('export PATH=/mnt/software/unstowable/anaconda/bin/:$PATH\n')
        fh.write('# snakemake:\n')
        # py3k env defined in /mnt/software/unstowable/anaconda and includes snakemake
        fh.write('source activate py3k;\n')
        #fh.write('cd {};\n'.format(os.path.abspath(args.outdir)))
        fh.write('cd $(dirname $0);\n'.format(os.path.abspath(args.outdir)))
        fh.write('# qsub for snakemake itself\n')
        fh.write('qsub="qsub -pe OpenMP 1 -l mem_free=4G -l h_rt=48:00:00 {} -j y -V -b y -cwd";\n'.format(mail_option))
        fh.write('# -j in cluster mode is the maximum number of spawned jobs\n')
        fh.write('$qsub -N vipr2.{} -o {}/snakemake.qsub.log'.format(conf['SAMPLENAME'], LOG_REL_DIR))
        qsub_per_task = "qsub -pe OpenMP {threads} -l mem_free=8G -l h_rt=24:00:00 -j y -V -b y -cwd"
        # FIXME max runtime and mem should be defined per target in SNAKEMAKE_FILE
        qsub_per_task += " -e {} -o {}".format(LOG_REL_DIR, LOG_REL_DIR)

        fh.write(' \'snakemake -j 8 -c "{}" -s {} --configfile {} --printshellcmds\';\n'.format(
            qsub_per_task, SNAKEMAKE_FILE, CONFIG_FILE))

    cmd = ['bash', snakemake_cluster_wrapper]
    if args.no_run:
        print("Not actually submitting job")
        print("When ready run: {}".format(' '.join(cmd)))
    else:
        print("Running: {}".format(cmd))
        subprocess.check_output(cmd)


if __name__ == '__main__':
    main()
