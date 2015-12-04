import shutil
from snakemake.utils import report

# simulate a bash login shell, see https://bitbucket.org/johanneskoester/snakemake/wiki/FAQ
shell.executable("/bin/bash")
# "unofficial bash strict mode" http://www.redsymbol.net/articles/unofficial-bash-strict-mode/
#shell.prefix("source ~/.bashrc; set -euo pipefail;")
shell.prefix("set -euo pipefail;")


from itertools import groupby

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
	input: 'assembly/assembly.fa'
	message: 'This is the end. My only friend, the end'
	output: 'COMPLETE'
	shell:  'touch {output}'

# trick is that SAMPLES was sorted in wrapper and can there safely be
# concatenated here. FIXME could we do the same here?

# famas will concatenate and write stats required for downsamling
rule concat_fastq:
	input:  fqs1=expand('{sample}R1.fastq.gz', sample=config['SAMPLES']),
			fqs2=expand('{sample}R2.fastq.gz', sample=config['SAMPLES'])
	output: fq1='R1.fastq.gz', fq2='R2.fastq.gz'
	log: "concat_fastq.log"
	message: "Concatenating FastQ files"
	benchmark: 'benchmark/concat_fastq.json'
	shell:  '{config[FAMAS]} -i <(zcat {input.fqs1}) -j <(zcat {input.fqs2}) -o {output.fq1} -p {output.fq2} 2>{log}'

rule downsample:
	input: fq1=rules.concat_fastq.output.fq1,
		   fq2=rules.concat_fastq.output.fq2,
		   reffa=config['REFFA']
	params: coverage='1000'
	benchmark: 'benchmark/downsample.json'
	log: 'downsample.log'
	message: "Downsampling"
	output: fq1='R1_1kcov.fastq.gz',
			fq2='R2_1kcov.fastq.gz'
	run:
		import sys
		with open(rules.concat_fastq.log[0]) as fh:
			# NOTE: assuming PE
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
	output: contigs='assembly/contigs.fasta'
	threads: 8
	log: 'assemble.log'
	benchmark: 'benchmark/assemble.json'
	shell:	"""
		#export MODULEPATH=""
		module use /mnt/software/modulefiles
		module load KMC smalt mummer samtools/0.0.19 python/2.7
		export CONDA_OLD_PS1="" PS1="" 
		source activate iva
		test -d ./assembly && rm -rf ./assembly
		iva -f {input.fq1} -r {input.fq2} -t {threads} ./assembly
		"""

rule join_contigs:
	input: contigs=rules.assemble.output.contigs, reffa=config['REFFA']
	output: assembly='assembly/assembly.fa'
	message: 'Joining contigs'
	params: samplename=config['SAMPLE_NAME'],
	log: 'join_contigs.log'
	benchmark: 'benchmark/join_contigs.json'
	shell: """
		simple_contig_joiner.py -c {input.contigs} -r {input.reffa} -o - | sed -e 's,^>joined,>{params.samplename}-assembly,' > {output}
		"""
	# NOTE: output might be empty of reference is too divergent
