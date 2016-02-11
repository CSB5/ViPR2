ViPR2: This is a remake of the original
[ViPR pipeline](https://github.com/CSB5/vipr) developed within the
[Genome Institute of Singapore](http://www.a-star.edu.sg/gis).

## How it works

ViPR2
- assembles your (downsampled) viral amplicons sequences using [IVA](http://www.ncbi.nlm.nih.gov/pubmed/25725497),
- orients assembled contigs fills gaps with the reference using [simple-contig-joiner](https://github.com/andreas-wilm/simple-contig-joiner),
- maps all reads with [BWA-MEM](http://arxiv.org/abs/1303.3997) against the given reference and the assembly
- determines primer positions to be ignored in next step
- and calls low frequency variants (SNVs and Indels) in both with [LoFreq](http://www.ncbi.nlm.nih.gov/pubmed/23066108)

The main script is called `vipr2.py`. Call it with `--help` to get some basic usage information.

After a successfull run you will find the main results in a folder
called results. The file `report.html` contains a description of all
relevant output files.


## Setup

ViPR2 is build with
[Snakemake](http://www.ncbi.nlm.nih.gov/pubmed/22908215) andq depends
on a multitude of software packages. The setup is heavily tuned for
our in-house settings. If you want to replicate it, start by changing
the CONF variable in `vipr2.py`, which lists expected binaries
etc. You will also have to change the way the jobs are submitted once
the config files and snakemake file have been created



