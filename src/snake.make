from snakemake.utils import report

shell.executable("/bin/bash")
# "unofficial bash strict mode" http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
shell.prefix("set -euo pipefail;")


from itertools import groupby


RESULT_DIR = "results"
RESULT_DIR_REFFA = os.path.join(RESULT_DIR, "reference-based")
RESULT_DIR_ASSEMBLY = os.path.join(RESULT_DIR, "assembly-based")
LOG_DIR = "logs"


localrules: vcf2csv, map_rate


def fasta_iter(fasta_name):
    """
    given a fasta file. yield tuples of header, sequence

    by Brent Pedersen: https://www.biostars.org/p/710/
    """
    fh = open(fasta_name)
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
        # drop the ">"
        header = next(header)[1:].strip()
        # join all sequence lines to one.
        seq = "".join(s.strip() for s in next(faiter))
        yield header, seq


rule final:
    input:
        os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.vcf.gz".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.csv".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.maprate.txt".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.covplot.pdf".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_REFFA, "{}.indelprep.vcf.gz".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_REFFA, "{}.indelprep.csv".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_REFFA, "{}.indelprep.maprate.txt".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR_REFFA, "{}.indelprep.covplot.pdf".format(config['SAMPLENAME'])),
        os.path.join(RESULT_DIR, "report.html")
    message: 'This is the end. My only friend, the end'


# famas will concatenate and write stats required for downsamling
rule concat_fastq:
    input:  fqs1=expand('{sample}R1.fastq.gz', sample=config['SAMPLES']),
            fqs2=expand('{sample}R2.fastq.gz', sample=config['SAMPLES'])
    output: fq1=temp('R1.fastq.gz'), fq2=temp('R2.fastq.gz')
    log: os.path.join(LOG_DIR, "concat_fastq.log")# needed for reading length and number of seqs later
    message: "Concatenating FastQ files"
    benchmark: 'benchmark/concat_fastq.json'
    shell:  '{config[FAMAS]} -i <(zcat {input.fqs1}) -j <(zcat {input.fqs2}) -o {output.fq1} -p {output.fq2} 2>{log}'


rule downsample:
    input: fq1=rules.concat_fastq.output.fq1,
           fq2=rules.concat_fastq.output.fq2,
           reffa=config['REFFA']
    params: coverage='1000'
    benchmark: 'benchmark/downsample.json'
    log: os.path.join(LOG_DIR, 'downsample.log')
    message: "Downsampling for assembly"
    output: fq1=temp('R1_1kcov.fastq.gz'),
            fq2=temp('R2_1kcov.fastq.gz')
    run:
        import sys
        with open(rules.concat_fastq.log[0]) as fh:
            #NOTE: assuming PE
            l = next(fh).strip()
            assert ' pairs in' in l
            l = next(fh).strip()
            assert ' pairs out'  in l
            num_pairs_out = int(l.split()[-1])
            l = next(fh).strip()
            assert 'length' in l
            avg_len = float(l.split()[-1])
        ref_seq = dict((x[0].split()[0], x[1])
                for x in fasta_iter(input.reffa))
        assert len(ref_seq) == 1, ("Only one reference sequence supported")
        ref_seq_len = len(list(ref_seq.values())[0])
        #sys.stderr.write("num_pairs_out={} avg_len={} ref_seq_len={}\n".format(num_pairs_out, avg_len, ref_seq_len))
        cov = num_pairs_out*2*avg_len / float(ref_seq_len)
        sample_every = int(cov / float(params.coverage))
        shell("{config[FAMAS]} --no-filter -i {input.fq1} -j {input.fq2} -o {output.fq1} -p {output.fq2} -s {sample_every} 2>{log}")


rule assemble:
    input: fq1=rules.downsample.output.fq1, fq2=rules.downsample.output.fq2
    output: contigs=os.path.join(RESULT_DIR_ASSEMBLY, 'assembly', 'contigs.fasta')
    params: outdir=os.path.join(RESULT_DIR_ASSEMBLY, 'assembly')
    threads: 8
    message: "Assembling downsampled reads with IVA"
    benchmark: 'benchmark/assembly.json'
    shell:  """
            {config[IVA]} {input.fq1} {input.fq2} {params.outdir} {threads}
            """

rule join_contigs:
    input: contigs=rules.assemble.output.contigs, reffa=config['REFFA']
    output: assembly=os.path.join(RESULT_DIR_ASSEMBLY, 'assembly', '{}-assembly.fa'.format(config['SAMPLENAME']))
    message: 'Joining contigs'
    params: samplename=config['SAMPLENAME']
    benchmark: 'benchmark/join_contigs.json'
    shell: """
        export PATH={config[MUMMERDIR]}:$PATH;
        {config[SIMPLE_CONTIG_JOINER]} -c {input.contigs} -r {input.reffa} -o - | sed -e 's,^>joined,>{params.samplename}-assembly,' > {output}
        """


rule bam_index:
    input:
        "{prefix}.bam"
    output:
        "{prefix}.bam.bai"
    #message: 'Creating BAM index for {input}'
    shell:
        "{config[SAMTOOLS]} index {input};"


# https://groups.google.com/forum/#!searchin/snakemake/lambda$20/snakemake/BidG8CEUW24/tsklilCYHJcJ
MAP_CMD = """
        {config[BWA]} index {input.reffa};# FIXME should be generic external rule
        rgid=$(echo {input.fq1} {input.fq2} | md5sum | cut -d " " -f1)
        rgpu=${{rgid}}.PU
        {config[BWA]} mem -M -t {threads} {input.reffa} {input.fq1} {input.fq2} -R "@RG\\tID:${{rgid}}\\tPL:illumina\\tPU:${{rgpu}}\\tSM:{config[SAMPLENAME]}" | \
            {config[SAMTOOLS]} sort -@ {threads} -o {output.bam} -
        """

rule map_reads_reffa:
    input:
        reffa=config['REFFA'],
        fq1=rules.concat_fastq.output.fq1,
        fq2=rules.concat_fastq.output.fq2
    output:
        bam=os.path.join(RESULT_DIR_REFFA, "{}.bam".format(config['SAMPLENAME']))
    message:'Mapping reads against {input.reffa}'
    threads: 8
    benchmark: 'benchmark/map_to_assembly.json'
    shell: 
        MAP_CMD

rule map_reads_assembly:
    input:
        reffa=rules.join_contigs.output.assembly,
        fq1=rules.concat_fastq.output.fq1,
        fq2=rules.concat_fastq.output.fq2
    output:
        bam=os.path.join(RESULT_DIR_ASSEMBLY, "{}.bam".format(config['SAMPLENAME']))
    message:'Mapping reads against {input.reffa}'
    threads: 8
    benchmark: 'benchmark/map_to_assembly.json'
    shell: 
        MAP_CMD

rule coverage_plot:
    input:
        bam='{prefix}.bam'
    output:
        covplot='{prefix}.covplot.pdf',
        covlog='{prefix}.covplot.txt'
    message:
        "Creating coverage plot"
    shell:
        """
        # for python2.7 with matplotlib
        export PATH=/mnt/software/unstowable/anaconda/bin/:$PATH;
        # for genomeCoverageBed
        export PATH=/mnt/software/stow/bedtools2-2.25.0/bin/:$PATH;
        # samtools
		export PATH=$(dirname {config[SAMTOOLS]}):$PATH;
        {config[COVERAGE_PLOT]} --force --log {output.covlog} -o {output.covplot} -b {input.bam};
        """


# taken from sg10k
rule samtools_fasta_index:
    input:
        "{prefix}.{suffix}"
    output:
        "{prefix}.{suffix,(fasta|fa)}.fai"
    #message: "Generating fasta index for {input}"
    shell:
        "{config[SAMTOOLS]} faidx {input};"


# taken from sg10k
rule map_rate:
    input:
        bam="{prefix}.idxstats.txt",
    output:
        "{prefix}.maprate.txt"
    #message: 'Computing mapping rate for {prefix}'
    shell:
        "cat {input} | awk '{{a+=$3; u+=$4}} END {{printf \"%.2f\\n\", a/(a+u)}}' > {output}"


# taken from sg10k
rule bam_idxstats:
    input:
        bam="{prefix}.bam",
        bai="{prefix}.bam.bai"
    output:
        "{prefix}.idxstats.txt"
    message: "Computing BAM stats for {input.bam}"
    shell:
        "{config[SAMTOOLS]} idxstats {input.bam} > {output};"


PRIMER_POS_CMD = """
        export PATH={config[MUMMERDIR]}:$PATH;
        # I'd love to use some other tools e.g. EMBOSS' primersearch but they are all equally unuseable
        {config[PRIMER_POS_FROM_SEQ]} -p {input.primer_fa} -r {input.reffa} --force -o {output.primerpos};
        seqname=$(cat {input.reffai} | cut -f 1)
        seqlen=$(cat {input.reffai} | cut -f 2)
        {config[PRIMER_POS_TO_BED]} --force -i {output.primerpos} --seqname $seqname --seqlen $seqlen -p {config[PRIMER_LEN]} -o {output.primersexclbed};
        """

rule determine_primer_pos_assembly:
    input:
        primer_fa=config['PRIMER_FILE'],
        reffa=rules.join_contigs.output.assembly,
        reffai=rules.join_contigs.output.assembly + ".fai"
    output:
        primerpos=os.path.join(RESULT_DIR_ASSEMBLY, "{}.primers.pos".format(config['SAMPLENAME'])),
        primersexclbed=os.path.join(RESULT_DIR_ASSEMBLY, "{}.primers.excl.bed".format(config['SAMPLENAME']))
    message:
        "Determining primer positions for {input.reffa}"
    shell:
        PRIMER_POS_CMD


rule determine_primer_pos_reffa:
    input:
        primer_fa=config['PRIMER_FILE'],
        reffa=config['REFFA'],
        reffai=config['REFFA'] + ".fai"
    output:
        primerpos=os.path.join(RESULT_DIR_REFFA, "{}.primers.pos".format(config['SAMPLENAME'])),
        primersexclbed=os.path.join(RESULT_DIR_REFFA, "{}.primers.excl.bed".format(config['SAMPLENAME']))
    message:
        "Determining primer positions for {input.reffa}"
    shell:
        PRIMER_POS_CMD


rule indel_calling_prep:
    """NOTE: could use lofreq viterbi instead of bamleftalign: that would fix alignment problems and do leftalignment (requires sorting), but multithreading option is missing  
    """
    input:
        bam="{prefix}.bam",
        reffa=lambda wildcards: config['REFFA'] if RESULT_DIR_REFFA in wildcards.prefix else rules.join_contigs.output.assembly
    output:
        bam="{prefix}.indelprep.bam"
    message:
        "Preparing BAM for indel calling"
    shell:
        "{config[LOFREQ]} indelqual --dindel -f {input.reffa} -o - {input.bam} | {config[BAMLEFTALIGN]} -c -f {input.reffa} > {output.bam}"


rule call_variants:
    input:
        bam='{prefix}.indelprep.bam',
        bai='{prefix}.indelprep.bam.bai',
        reffa=lambda wildcards: config['REFFA'] if RESULT_DIR_REFFA in wildcards.prefix else rules.join_contigs.output.assembly,
        reffai=lambda wildcards: config['REFFA'] + ".fai" if RESULT_DIR_REFFA in wildcards.prefix else rules.join_contigs.output.assembly + ".fai",
        nonprimer_regions=lambda wildcards: rules.determine_primer_pos_reffa.output.primersexclbed  if RESULT_DIR_REFFA in wildcards.prefix else rules.determine_primer_pos_assembly.output.primersexclbed
    threads:
        8
    output:
        vcf='{prefix}.indelprep.vcf.gz'
    shell:
        """export PATH=/mnt/software/bin/:$PATH
        {config[LOFREQ]} call-parallel --pp-threads {threads} -f {input.reffa} -l {input.nonprimer_regions} --call-indels -o {output.vcf} {input.bam}
        """
    message:
        "Calling variants with LoFreq"


rule vcf2csv:
    input:
        "{prefix}.vcf.gz"
    output:
        "{prefix}.csv"
    shell:  
        """# for python2.7 with pyvcf
        export PATH=/mnt/software/unstowable/anaconda/bin/:$PATH;
        {config[VCF2CSV]} {input} {output}
        """
    message:
        "Converting vcf to csv"


rule report:
    # ugly. rules.call_variants.output.vcf doesnt work (can't find 'sample')
    input: assembly_reffa=rules.join_contigs.output.assembly,
           assembly_vcf=os.path.join(RESULT_DIR_ASSEMBLY, '{}.indelprep.vcf.gz'.format(config['SAMPLENAME'])),
           assembly_csv=os.path.join(RESULT_DIR_ASSEMBLY, '{}.indelprep.csv'.format(config['SAMPLENAME'])),
           assembly_covplot=os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.covplot.pdf".format(config['SAMPLENAME'])),
           assembly_maprate=os.path.join(RESULT_DIR_ASSEMBLY, "{}.indelprep.maprate.txt".format(config['SAMPLENAME'])),
           reference_vcf=os.path.join(RESULT_DIR_REFFA, '{}.indelprep.vcf.gz'.format(config['SAMPLENAME'])),
           reference_csv=os.path.join(RESULT_DIR_REFFA, '{}.indelprep.csv'.format(config['SAMPLENAME'])),
           reference_covplot=os.path.join(RESULT_DIR_REFFA, "{}.indelprep.covplot.pdf".format(config['SAMPLENAME'])),
           reference_maprate=os.path.join(RESULT_DIR_REFFA, "{}.indelprep.maprate.txt".format(config['SAMPLENAME'])),
    output: html=os.path.join(RESULT_DIR, "report.html")
    run: report("""
          ===================
          ViPR2 Report
          ===================

          ViPR2 will assemble your (downsampled) viral amplicon sequencing reads using 
          IVA (http://www.ncbi.nlm.nih.gov/pubmed/25725497). Contigs are oriented and gaps
          are filled with the provided reference sequence. Reads are mapped against this
          assembly and the input reference with BWA-MEM (http://arxiv.org/abs/1303.3997). 
          Low-frequency SNVs and Indels are then called (ignoring determined primer positions) 
          with LoFreq (http://www.ncbi.nlm.nih.gov/pubmed/23066108).

          Output files
          ------------

          - Assembled sequence: {input.assembly_reffa}
          - Variants wrt. assembly: {input.assembly_vcf} or for Excel import {input.assembly_csv}
          - Coverage plot for assembly mapping: {input.assembly_covplot}
          - Mapping success: {input.assembly_maprate}

          - Variants wrt. reference: {input.reference_vcf} or for Excel import {input.reference_csv}
          - Coverage plot for reference mapping: {input.reference_covplot}
          - Mapping success: {input.reference_maprate}

          See conf.json in the main directory for used programs and settings.
          """, output.html, metadata="Andreas WILM", **input)

