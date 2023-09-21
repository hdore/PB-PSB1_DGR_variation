# Scripts used in the PB-PSB1 DGR article

This github repository holds the scripts used to perform the analyses described in the article [insert title and DOI and link].

## Python scripts

### Calculation of the number of potential variants generated by DGR in a target

The python script [count_AA_possibilities_in_DGR_regions_v3.py](./python_scripts/count_AA_possibilities_in_DGR_regions_v3.py) calculates the number of protein variants that can be generated by one or multiple targets, working at the codon level.

The script takes a fasta file with a single sequence (it does not support multi-fasta files, e.g. genomes with multiple contigs) and a table of DGR coordinates with target, VR and TR coordinates and TR sequence, and outputs the number of possible amino-acids for each putatively variable position in each VR of each target. It takes into account the cases where multiple variable positions are present within a codon. The codon sequence is taken from the TR, but the VR coordinates are used to determine the coding frame. 

The script can either extract the TR sequence from the TR coordinates (default), or from the TR sequence indicated in the input file (with the flag `--tr_from_file`) to accomodate potential VR-TR indels (in our case a 3bp deletion in the VR).

#### Input files

For the input fasta file, only single-sequence fasta is supported. If your genome has multiple contigs, please use a fasta file with only the contig of interest.

The input "regions file" with regions coordinates should be formatted as in the [example_regions_file_for_count_AA_possibilities](./example_input_files/example_regions_file_for_count_AA_possibilities.txt). If you include a header in this file, it must be commented with a # (see example file). Note that to be accurate, the coordinates should be 1-based (as found in a genbank file), i.e. the start coordinate should correspond to the first base, and the end coordinate to the last base. 
**Important: if the Target is on the minus strand, the start coordinate must be greater than the end coordinate (this is how the script recognizes the strand, allowing it to extract the sequence correctly). The same applies to VR and TR coordinates.**

#### Output

The tab-delimited output file has a line for each DGR-targeted codon and 14 columns:

- CodonID: a unique identifier for the codon.

- DGR: the name of the DGR locus (from the input table of coordinates)

- TargetID: the name of the target (from the input table of coordinates)

- VRnumber: the name of the VR within a locus (from the input table of coordinates)

- VRid: a unique identifier for the VR (from the input table of coordinates)

- VR_start:  the start position of the VR (1-based)

- VR_end: the end position of the VR (1-based, inclusive)

- A_pos_in_target_1based: the position corresponding to a TR adenine within the target gene

- codon_start_1based: the coordinate of the first base of the codon

- ref_codon: the codon in the reference sequence of the TR

- place_in_codon: the place in the codon where adenines are found in the TR (1 and/or 2 and/or 3, comma separated)

- nb_possible_AA: the number of possible amino-acid encoded by this codon when adenines are replaced by any nucleotide

- nb_possible_stop: the number of possible stop codons generated when adenines are replaced by any nucleotide

- VR_nb_possible_combinations: the number of possible combinations of amino-acids at the VR level (this will be the same value for all codons of a given VR)

#### Script usage

To see the script usage, run `python count_AA_possibilities_in_DGR_regions_v3.py -h`:

[paste usage]

```
usage: count_AA_possibilities_in_DGR_regions_v3.py [-h] [-r REGIONS_FILE]
                                                   [-o OUTPUT_FILE]
                                                   [-i --input_file INPUT_FILE]
                                                   [--tr_from_file] [-d]

Takes a fasta file with a single sequence and a table of DGR coordinates with
target, VR and TR coordinates and TR sequence, and outputs the number of
possible amino-acids for each putatively variable position in each VR of each
target. V1 implements a way to take into account multiple variable positions
within a codon. v2 considers codon sequence from the TR, not the target (but
still using target sequence to know where codons start and end). v3 adds the
possibility to use TR sequence from the input file to accomodate potential
mismatches (in our case a 3bp deletion in the VR), and adds codon sequence to
output file.

optional arguments:
  -h, --help            show this help message and exit
  -r REGIONS_FILE, --regions_file REGIONS_FILE
                        Name or path to input file with the VR regions
                        coordinates. Must be 'DGR Target_ID VR_number VR_ID
                        Target start Target end VR start VR end VRsize TR
                        start TR end TR seq'.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Name or path for output tab-delimited file.
  -i --input_file INPUT_FILE
                        Name or path of the fasta file of the genome (for now
                        only single-sequence fasta files are accepted).
  --tr_from_file        Flag. If set, the program will use the TR sequence
                        provided in the input regions file. Otherwise, it
                        extracts the TR sequence from the fasta file based on
                        the TR coordinates. By default the flag is not set and
                        the sequence is extracted based on TR coordinates.
  -d, --debug_mode      Enable debug mode (very verbose).
```

#### Dependencies

This script makes use of [Biopython](https://biopython.org). Biopython must be present on the system, and can be installed with conda: `conda install biopython` 

### Filter bam file to get read pairs mapping on DGR variable regions

The python script [filter_bam_by_mate_location.py](./python_scripts/filter_bam_by_mate_location.py) allows to filter a bam file to keep only pairs of reads in which one of the reads aligns to a variable region, and the other read of the pair is within a certain genomic distance set by the user.

#### Input files

The script takes as input a list of sorted and indexed bam files (1 per line), a file with start and end coordinates of regions of interest (see [example_regions_file_for_filter_bam_and_calc_pi](./example_input_files/example_regions_file_for_filter_bam_and_calc_pi.txt)). If you include a header in this file, it must be commented with a #. 

The script was only tested with bam files containing a single reference (i.e. with reads mapped to a single-sequence fasta file). It is very likely that it won't work with bam files containing multiple references. If you have a genome with multiple contigs, choose only the contig containing the DGR locus or loci before doing the mapping, or extract this contig from the bam file before running the script.

#### Output

The script outputs new bam files (as many as in the input list), with only reads mapping within the region of interest whose mate is within the defined window (1000 bp by default).

#### Script usage

To see all arguments, run `python filter_bam_by_mate_location.py -h`

```
usage: filter_bam_by_mate_location.py [-h] [-b BAM_LIST] [-r REGIONS_FILE]
                                      [-s SUFFIX]
                                      [--mq_threshold MQ_THRESHOLD]
                                      [--window WINDOW] [-d]

Takes as input a list of bam files (1 per line), start and end coordinates of
regions of interest and a window size, and creates new bam files with only
reads mapping within the region of interest whose mate is within the defined
window.

optional arguments:
  -h, --help            show this help message and exit
  -b BAM_LIST, --bam_list BAM_LIST
                        Name or path to input list of bam files. Must be 1
                        file/line, with full path.
  -r REGIONS_FILE, --regions_file REGIONS_FILE
                        Name or path to input file with the VR regions
                        coordinates.
  -s SUFFIX, --suffix SUFFIX
                        Suffix for output subsetted bam file. The '.bam' of
                        the input bam will be replaced by 'suffix.bam'.
                        Default is 'subsetted'
  --mq_threshold MQ_THRESHOLD
                        Int. Minimum mapping quality for a read to be
                        considered. Default is 0.
  --window WINDOW       Int. Window in which the mate is looked for. The mate
                        will be looked for in the sequence defined by
                        [region_start - window, region_end + window]. Default
                        is 1000.
  -d, --debug_mode      Enable debug mode (very verbose).
```

#### Dependencies

This script makes use of [Biopython](https://biopython.org) and [pysam](https://pypi.org/project/pysam/). Both must be present on the system, and can be installed with conda: `conda install biopython ; conda install -c bioconda pysam`.

### Calculation of per-position nucleotide diversity and proportion of non-reference alleles from a metagenome

The python script [calc_pi_from_pileup_v1.py](./python_scripts/calc_pi_from_pileup_v1.py) allows to calculate the nucleotide diversity ($\pi$) and the proportion of non-reference alleles from a bam file pileup, at positions defined by the user.

The nucleotide diversity (pi) can be seen as the probability for 2 reads to present the same nucleotide at a given position. It is calculated as $\pi=1-(a^2+c^2+t^2+g^2)$ with a, t, c, g the frequency of each nucleotide at the position of interest (see also the Methods in [Soil bacterial populations are shaped by recombination and gene-specific selection across a grassland meadow | The ISME Journal](https://doi.org/10.1038/s41396-020-0655-x)).

The proportion of non-reference alleles (nonref_ratio) is the proportion of reads at a given position that show a nucleotide differing from the reference nucleotide.

Both metrics are calculated for each individual position.

#### Input files

The script works on a pileup generated from a (sorted and indexed) bam file. It was only tested with bam files containing a single reference (i.e. with reads mapped to a single-sequence fasta file). It is very likely that it won't work with bam files containing multiple references. If you have a genome with multiple contigs, choose only the contig containing the DGR locus or loci before doing the mapping, or extract this contig from the bam file before running the script.

To generate a pileup file from a bam file, you should use `samtools mpileup` with options `-a --fasta-ref your-reference.fasta`. For example:

```bash
samtools mpileup -a --fasta-ref /your_path/to_reference_genome/reference_genome.fasta your_input_sorted_indexed_bam.sorted.bam > your_mpileup.tab
#These two options are needed: --fasta-ref to have the reference allele and -a to output positions with no read
```

The script also needs a "regions file" containing the coordinates of the regions of interest (e.g. DGR variable repeats). The coordinates file must be a tab-delimited file with the following columns (see [example_regions_file_for_filter_bam_and_calc_pi](./example_input_files/example_regions_file_for_filter_bam_and_calc_pi.txt)): 
`#Target_ID    VR_number    VR_ID    Target start    Target end    VR start    VR end`

If you include a header in this file, it must be commented with a # (see example file). 
**Important: the VR start coordinate must always be lower than the VR end coordinate.** (I know this can be confusing between my different scripts, and might fix it in the future).

#### Output

(Note that this script was designed for the study of an organism with many DGR loci. For this reason, it gives a condensed output that is not the easiest to parse. I you have a single DGR locus, consider using the script `nt_prop_from_pileup_v2.py` below, which also outputs the nucleotide diversity and proportion of non-reference alleles, with one line per genomic position.)

Here the output is a tab-delimited file, with the following colums:

`VRid    target_id    VR_start     VR_end    mean_pi    mean_nonref_ratio    pi_list    nonref_ratio_list    n_and_del_count_list`

- VRid: a unique identifier for the VR (from the input "regions file")

- targetID: the name of the target (from the input "regions file")

- VR_start: the start coordinate for this VR (from the input "regions file")

- VR_end: the end coordinate for this VR (from the input "regions file")

- mean_pi: the mean of pi values for all positions in this VR passing the coverage threshold

- mean_nonref_ratio: the mean of nonref_ratio values for all positions in this VR passing the coverage threshold

- pi_list: a comma-separated list of values of the nucleotide diversity (pi) for each position of the VR, with -1 if the position did not pass the coverage threshold

- nonref_ratio_list: a comma-separated list of values of the nonref_ratio for each position of the VR, with -1 if the position did not pass the coverage threshold

- n_and_del_count_list: a comma-separated list of the number of N or deletions found at each position of the VR

#### Script usage

To see all options and arguments, run `python calc_pi_from_pileup_v1.py -h`:

```
usage: calc_pi_from_pileup_v1.py [-h] [-r REGIONS_FILE] [-o OUTPUT_FILE]
                                 [-i --input_file INPUT_FILE]
                                 [-c --cov_cutoff COV_CUTOFF] [-d]

This script allows to calculate the nucleotide diversity (pi) and the
proportion of non-reference alleles from a bam file pileup, at positions
defined by the user.

optional arguments:
  -h, --help            show this help message and exit
  -r REGIONS_FILE, --regions_file REGIONS_FILE
                        Name or path to input file with the VR regions
                        coordinates.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Name or path for output tab-delimited file.
  -i --input_file INPUT_FILE
                        Name or path of input_file from samtools mpileup with
                        -a and --fasta-ref options.
  -c --cov_cutoff COV_CUTOFF
                        Int. The coverage cutoff at a position to calculate
                        statistics.
  -d, --debug_mode      Enable debug mode (very verbose).
```

#### Dependencies

This script makes use of [Biopython](https://biopython.org/) and [numpy](https://numpy.org/). Both must be present on the system, and can be installed with conda: `conda install biopython numpy`.

### Extraction of the proportion of each nucleotide at VR positions from a metagenome

The python script [nt_prop_from_pileup_v2.py](./python_scripts/nt_prop_from_pileup_v2.py) allows to report the frequency of each nucleotide (a, t, c, g), the nucleotide diversity ($\pi$) and the proportion of non-reference alleles, from a bam file pileup, at positions defined by the user.

The nucleotide diversity (pi) can be seen as the probability for 2 reads to present the same nucleotide at a given position. It is calculated as \pi=1-(a^2+c^2+t^2+g^2) with a, t, c, g the frequency of each nucleotide at the position of interest (see also the Methods in [Soil bacterial populations are shaped by recombination and gene-specific selection across a grassland meadow | The ISME Journal](https://doi.org/10.1038/s41396-020-0655-x)).

The proportion of non-reference alleles (nonref_ratio) is the proportion of reads at a given position that show a nucleotide differing from the reference nucleotide.

Both metrics are calculated for each individual position.

#### Input files

The script works on a pileup generated from a (sorted and indexed) bam file. It was only tested with bam files containing a single reference (i.e. with reads mapped to a single-sequence fasta file). It is very likely that it won't work with bam files containing multiple references. If you have a genome with multiple contigs, choose only the contig containing the DGR locus or loci before doing the mapping, or extract this contig from the bam file before running the script.

To generate a pileup file from a bam file, you should use `samtools mpileup` with options `-a --fasta-ref your-reference.fasta`. For example:

```bash
samtools mpileup -a --fasta-ref /your_path/to_reference_genome/reference_genome.fasta your_input_sorted_indexed_bam.sorted.bam > your_mpileup.tab
#These two options are needed: --fasta-ref to have the reference allele and -a to output positions with no read
```

The script also needs a "regions file" containing the coordinates of the regions of interest (e.g. DGR variable repeats). The coordinates file must be a tab-delimited file with the following columns (see [example_regions_file_for_nt_prop](./example_input_files/example_regions_file_for_nt_prop.txt)): `#Target_ID    VR_number    VR_ID    Target start    Target end    VR start    VR end    TR start    TR end`

If you include a header in this file, it must be commented with a # (see example file). 
**Important: if the Target is on the minus strand, the start coordinate must be greater than the end coordinate (this is how the script recognizes the strand, allowing it to extract the sequence correctly). The same applies to VR and TR coordinates.**

#### Output

The output is a single tab-delimited file with the following columns:

`VRid    target_id    VR_start    VR_end    VR_strand    position    ref_allele    TR_nt    nt_value    nt_count    pi    nonref_ratio`

This is a "long" format designed to be easily plotted with ggplot in R, so for each position there are 5 lines (for A count, C count, T count, G count and the count of Ns or indels. The corresponding nucleotide is indicated in the column nt_value). The pi and nonref_ratio values will be the same for these 5 lines,  since they are calculated by position.

- VRid: a unique identifier for the VR (from the input "regions file")

- target_id: the name of the target (from the input "regions file")

- VR_start: the start coordinate for this VR (from the input "regions file")

- VR_end: the end coordinate for this VR (from the input "regions file")

- position: the genomic position for which information is reported on this line

- ref_allele: the nucleotide displayed by the reference at this position

- TR_nt: the nucleotide displayed by the TR in the reference genome, at the position that aligns to the reported position

- nt_value: the nucleotide for which the count is reported on this line in the column nt_count (A or C
   or T or G or "Ns + indels")

- nt_count: the read count (i.e the number of reads having the nucleotide indicated in the column nt_value)

- pi: the nucleotide diversity calculated for this position. Note that the same value is displayed for the 5 lines describing each position.

- nonref_ratio: the proportion of non-reference alleles for this position. Note that the same value is displayed for the 5 lines describing each position.

#### Script usage

To see all options and arguments, run `python nt_prop_from_pileup_v2.py -h`:

```
usage: nt_prop_from_pileup_v2.py [-h] [-r REGIONS_FILE] [-o OUTPUT_FILE]
                                 [-i INPUT_FILE] [-c COV_CUTOFF] [-d]

Takes a table file with the coordinates of the regions of interest, and a
pileup file obtained by running samtools mpileup on a bam file (with -a and
--fasta-ref options). Outputs a table of each nt proportion at each position
of each region. It also calculates nucleotide diversity (pi) and non-reference
allele proportion as a bonus. v1 modifies pi and non-reference allele
proportion to use the sum of A + T + C + G instead of the coverage calculated
by mpileup, which includes Ns. V2 brings a correction for when VR is on the
minus strand.

optional arguments:
  -h, --help            show this help message and exit
  -r REGIONS_FILE, --regions_file REGIONS_FILE
                        Name or path to input file with the VR and TR regions
                        coordinates.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Name or path for output tab-delimited file.
  -i INPUT_FILE, --input_file INPUT_FILE
                        Name or path of input_file from samtools mpileup with
                        -a and --fasta-ref options.
  -c COV_CUTOFF, --cov_cutoff COV_CUTOFF
                        Int. The coverage cutoff at a position to calculate
                        statistics. Default: 5.
  -d, --debug_mode      Enable debug mode (very verbose).
```

#### Dependencies:

This script makes use of [Biopython](https://biopython.org/) and [numpy](https://numpy.org/). Both must be present on the system, and can be installed with conda: `conda install biopython numpy`.

## 
