#!/bin/env python2

import sys
import csv

import vcf


def vcf2csv(vcffile, csvfile):
    """FIXME:add-doc
    """

    
                                
    vcffh = vcf.VCFReader(filename=vcffile)
    fieldnames = ['Seq', 'Pos', 'Ref.Base', 'Var.Base', 'Qual', 'Allele Freq', 'Type', 'Homopolymer Length', 'Depth', 'Depth: ref fw and rev, var fw and rev']
    if csvfile == "-":
        csvfile = '/dev/stdout'
    with open(csvfile, 'wb') as csvfh:
        csvw = csv.DictWriter(csvfh, fieldnames=fieldnames)#, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        csvw.writeheader()
        for var in vcffh:
            assert len(var.ALT) == 1
            fields = [var.CHROM, var.POS, var.REF, var.ALT[0], var.QUAL, var.INFO['AF']]

            if var.is_indel:
                fields.extend(["INDEL", var.INFO['HRUN']])
            else:
                fields.extend(["SNV", "None"])
            fields.append(var.INFO['DP'])
            fields.append(','.join(str(x) for x in var.INFO['DP4']))
            row = dict(zip(fieldnames, fields))
            csvw.writerow(row)

def main():
    """main function
    """

    assert len(sys.argv) == 3, ("Need one vcf file and one csv file as arguments")
    vcffile = sys.argv[1]
    csvfile = sys.argv[2]
    vcf2csv(vcffile, csvfile)


if __name__ == '__main__':
    main()
