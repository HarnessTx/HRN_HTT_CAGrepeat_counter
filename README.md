# HRN_CAG_repeat_counter
HTT CAG Repeat Counter is a tool developed by Harness Therapeutics Limited for calculating the number of CAG repeats in the human HTT gene using nanopore sequencing data from PCR-amplified gDNA. Designed for analyzing repeat expansion disorders, this script detects HTT CAG repeat number from filtered, mapped fastq file.

## Description

The HRN_CAG_repeat_counter detects CAG repeat number in the human HTT gene, for each read in the input. It takes as input nanopore sequencing data from PCR amplified DNA, where the PCR amplicon spans the HTT exon 1 CAG repeat tract, and the interrupting motif region. The interrupting motif is defined as (CAACAG)<sub>1</sub>(CCGCCA)<sub>1</sub>. We code the interrupting motif as IUPAC `CARCARCCRCCR`.

The sequence of HTT exon 1 in the region of interest can be shown as (CAG)<sub>n</sub>(CAACAG)<sub>1</sub>(CCGCCA)<sub>1</sub>(CCG)<sub>7</sub>(CCT)<sub>2</sub>, where n=19 if it is from the reference allele. 
HRN_CAG_repeat_counter filters the HTT mapped reads for those containing the user provided PCR primer sequences, searches for the presence of the interrupting motif. The position of the forward primer relative to the start of the CAG repeat tract and to define the beginning of the CAG repeat. The script calculates the base pair (bp) distance between this position and the interrupting motif. This distance is then divided by 3 to determine the number of CAG repeats.

## Requirements

- R (only tested with version 4.4.3)
- Required R packages: `Biostrings` (Pagès _et al._, 2024), `ShortRead` (Morgan _et al._, 2009), `optparse`, `openxlsx`

## Installation

Clone the repository:

```bash
git clone https://github.com/HarnessTx/HRN_HTT_CAGrepeat_counter.git
cd HRN_HTT_CAGrepeat_counter
```

## Usage

```bash
Rscript HTT_CAG_counter.R --mapped_fq <path-to-mapped-fq> --name <name> --primer_seq_file <path-to-primer-csv>
```

## Run test dataset

```bash
Rscript HTT_CAG_counter.R --mapped_fq example_data/HTT_mapped_reads.fq.gz --name testdataset --primer_seq_file example_data/PrimerSequence_Goold.2021.csv
```

## Input data

This tool takes mapped reads in fastq format as the input. We suggest the following workflow for pre-processing raw nanopore reads:
1. Raw reads with an average quality < 12 are removed using nanoq (v0.10.0; Steinig and Coin 2022) 
2. Filtered reads are mapped to a pseudo-genome consisting of HTT exon 1 with CAG<sub>n</sub>, where n = 19, 109, 127, 144 or any number of the estimated CAG size on the mutant allele, using minimap2 (v2.26-r1175; Li 2018). 
3. SAM files are converted to BAM files and split into mapped and unmapped FASTQ files using samtools (v1.19.2; Li _et al._, 2009).
```
samtools view -b -F 4 aligned.bam | samtools bam2fq - | gzip > HTT_mapped_reads.fq.gz
```
4. HRN_CAG_repeat_counter requires one input fastq file - containing reads that were mapped to HTT.

## Required arguments
-m MAPPED_FQ, --mapped_fq=MAPPED_FQ
    Path to the input FASTQ file containing reads mapping to the pseudo-genome containing human HTT gene, as described in item 2 under `Input data`.

-p PRIMER_SEQ_FILE, --primer_seq_file=PRIMER_SEQ_FILE
    Path to the input csv file containing the primer sequences used for PCR amplification of human HTT gene

## Optional arguments
-n NAME, --name=NAME
    Name to append to output files

-h, --help
    Show help message and exit

## Output files

### HTTcagEstimationTable.xlsx
An Excel (xlsx) file containing the positions of detected interrupting motifs in each read.

- **names** - ID of the nanopore sequencing read taken from the input fastq.
- **Str.INTmotif** - Start position of the interruping motif, within the read.
- **End.INTmotif** - End position of the interrupting motif, within the read.
- **Width.INTmotif** - Width of the interrupting motif.
- **InterruptMotif.cat** - Indicates if the interrupting motif was detected in the FR or RC orientation.
- **Str.5prime** - Start position of the 5' primer, within the read.
- **End.5prime** - End position of the 5' primer, within the read.
- **Width.5prime** - Width of the 5' primer sequence.
- **Str.3prime** - Start position of the 3' primer, within the read.
- **End.3prime** - End position of the 3' primer, within the read.
- **Width.3prime** - Width of the 3' primer sequence.
- **EitherPrimer.cat** - Indicates the FR or RC orientation of the read.
- **QualityControl** - Indicates whether the read was mapped to the HTT gene
- **Dist.prim2INTmotif** - Distance between the end of 5' primer to the start of the interrupting motif.
- **Num.Trinuc** - Number of trinucleotides from `Dist.prim2INTmotif` to infer the number of CAG from HTT gene.


## Reference
>Goold R. _et al._. (2021). FAN1 controls mismatch repair complex assembly via MLH1 retention to stabilize CAG repeat expansion in Huntington's disease. _Cell Rep._ 36(9):109649. https://doi.org/10.1016/j.celrep.2021.109649. 

>Li H. _et al._ (2009). 1000 Genome Project Data Processing Subgroup. The Sequence Alignment/Map format and SAMtools. _Bioinformatics_. 25(16):2078-9. https://doi.org/10.1093/bioinformatics/btp352.

>Li H. (2018). Minimap2: pairwise alignment for nucleotide sequences. _Bioinformatics_. 34(18):3094-3100. https://doi.org/10.1093/bioinformatics/bty191.

>Morgan M. _et al._ (2009) ShortRead: a bioconductor package for input, quality assessment and exploration of high-throughput sequence data. _Bioinformatics_. 25(19):2607-8. https://doi.org/10.1093/bioinformatics/btp450.

>Pagès H. _et al._ (2024). Biostrings: Efficient manipulation of biological strings. https://doi.org/10.18129/B9.bioc.Biostrings, _R package version 2.74.1_. https://bioconductor.org/packages/Biostrings.

>Steinig and Coin (2022). Nanoq: ultra-fast quality control for nanopore reads. _Journal of Open Source Software_. 7(69), 2991. https://doi.org/10.21105/joss.02991.


