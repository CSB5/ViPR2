#!/mnt/software/unstowable/anaconda/bin/python2
"""The script will take a BAM file and fw/rv primer positions and
remove any forward reads starting at forwards positions and reverse
reads ending at the reverse positions

Harcoded python interpreter so that pysam hack below works
"""

#--- standard library imports
#
import os
import sys
import logging
import argparse
import gzip

#--- third-party imports
#
sys.path.insert(0, "/mnt/software/unstowable/anaconda/envs/pysam-0.8.4/lib/python2.7/site-packages")
import pysam
assert [int(x) for x in pysam.__version__.split(".")] >= [0, 8, 0]

#--- project specific imports
#
import primer_pos as pp


__author__ = "Andreas Wilm"
__version__ = "0.2"
__email__ = "wilma@gis.a-star.edu.sg"
__copyright__ = "2012, 2013, 2014, 2015, 2016 Genome Institute of Singapore"
__license__ = "GPL2"
__status__ = "eternal alpha"


#global logger
# http://docs.python.org/library/logging.html
LOG = logging.getLogger("")
logging.basicConfig(level=logging.WARN,
                    format='%(levelname)s [%(asctime)s]: %(message)s')


def cmdline_parser():
    """
    creates an argparse instance
    """

    parser = argparse.ArgumentParser(__doc__)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Optional: be verbose")
    parser.add_argument("--force",
                        action="store_true",
                        dest="force_overwrite",
                        help="Optional: force overwriting of files")
    parser.add_argument("--debug",
                        action="store_true",
                        dest="debug",
                        help="Optional: enable debugging")
    parser.add_argument("-i", "--bam-in",
                        dest="fbam_in",
                        help="BAM input file (stdin not supported)")
    parser.add_argument("-o", "--bam-out",
                        dest="fbam_out",
                        help="BAM output file (stdout supported)")
    parser.add_argument("-p", "--primer",
                        dest="primer_pos_file",
                        help="Primer start position file Format:"
                        " 1-based-start-pos orientation-F-or-R.")
    return parser



def fragment_name(r):
    """
    >>> fragment_name("@EAS139:136:FC706VJ:2:2104:15343:197393 1:N:18:1")
    '@EAS139:136:FC706VJ:2:2104:15343:197393'

    >>> fragment_name("@HWUSI-EAS100R:6:73:941:1973#0/1")
    '@HWUSI-EAS100R:6:73:941:1973#0'
    """
    
    if r[-2] == "/":
        return r[:-2]
    else:
        return r.split()[0]
    

        
def find_primer_reads(bam_in, peaks_fw_start_pos, peaks_rv_end_pos, dups_fh):
    """Find primer reads
    """

    bam_fh = pysam.Samfile(bam_in, "rb")
    for r in bam_fh:
        if r.is_unmapped or r.is_secondary or r.is_supplementary:
            continue

        is_dup = False
        if r.is_reverse:
            if r.reference_end-1 in peaks_rv_end_pos:
                is_dup = True
        else:
            if r.reference_start in peaks_fw_start_pos:
                is_dup = True

        if is_dup:
            dup_frag_name = fragment_name(r.query_name)
            dups_fh.write(dup_frag_name + "\n")
                
    bam_fh.close()


def mark_primer(sam_in, sam_out, dups_fh):
    """FIXME:add-doc
    """

    dup_names = dict()
    num_marked_as_dups = 0

    for line in dups_fh:
        dup_names[line.rstrip()] = 1
            
    LOG.info("Names of {} duplicated fragments read from {}".format(
        len(dup_names), dups_fh.name))
    for r in sam_in:
        #LOG.critical("{}, {}...".format(fragment_name(r.query_name), dup_names.keys()[:10]))
        if dup_names.has_key(fragment_name(r.query_name)):
            #LOG.critical("R was {}".format(r))
            r.flag |= 0x400
            #LOG.critical("R is {}".format(r))
            num_marked_as_dups += 1
        sam_out.write(r)

    LOG.info("Marked {} reads as dups".format(num_marked_as_dups))


def main():
    """The main function
    """

    parser = cmdline_parser()
    args = parser.parse_args()

    if args.verbose:
        LOG.setLevel(logging.INFO)
    if args.debug:
        LOG.setLevel(logging.DEBUG)

    # file check
    #
    for (filename, descr, direction, mandatory) in [
            (args.fbam_in, 'BAM input file', 'in', True),
            (args.fbam_out, 'BAM output file', 'out', True),
            (args.primer_pos_file, 'Primer positions file', 'in', True)]:

        if not mandatory and not filename:
            continue

        if not filename:
            parser.error("%s argument missing." % descr)
            sys.exit(1)

        if filename == '-':
            continue

        if direction == 'in' and not os.path.exists(filename):
            LOG.fatal("file '{}' does not exist.\n".format(filename))
            sys.exit(1)

        if direction == 'out' and os.path.exists(filename) and not args.force_overwrite:
            LOG.fatal("Refusing to overwrite existing file '{}'.\n".format(filename))
            sys.exit(1)

    if args.fbam_in == "-":
        LOG.fatal("No streaming allow for input BAM (read twice)")
        sys.exit(1)

    # parse peaks and prep
    #
    peaks = pp.parse_primer_pos(open(args.primer_pos_file))
    # summarize fw-start and rv-end positions in list. actually don't
    # need anything else (i.e. the end positions)
    peaks_fw_start_pos = [p.pos for p in peaks if p.ori == 'F']
    peaks_rv_end_pos = [p.pos for p in peaks if p.ori == 'R']

    dupreads_out = args.fbam_out.replace(".bam", ".dups.gz")
    if os.path.exists(dupreads_out):
        LOG.warn("Reusing {}".format(dupreads_out))
    else:
        with gzip.open(dupreads_out, 'w') as dupreads_fh:
            find_primer_reads(args.fbam_in, peaks_fw_start_pos, peaks_rv_end_pos, dupreads_fh)

    with pysam.Samfile(args.fbam_in, "rb") as sam_in:
        with pysam.Samfile(args.fbam_out, "wb", template=sam_in) as sam_out:
            with gzip.open(dupreads_out) as dupreads_fh:
                mark_primer(sam_in, sam_out, dupreads_fh)


if __name__ == "__main__":
    main()
    LOG.info("Successful exit")
