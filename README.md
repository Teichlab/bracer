# BraCeR
BraCeR - reconstruction of B cell receptor sequences from single-cell RNA-seq data.

**IMPORTANT: Python dependencies are different from TraCeR. Use the requirements file (detailed [here](#setup)) to update them.**

## Contents ##
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Setup](#setup)
4. [Testing](#testing-bracer)
5. [Usage](#using-bracer)
	- [*Assemble*](#assemble-bcr-reconstruction)
    - [*Summarise*](#summarise-summary-and-clonotype-networks)
6. [Docker image](#docker-image)


## Introduction
This tool reconstructs the sequences of rearranged and expressed B cell receptor genes from single-cell RNA-seq data. It then uses the BCR sequences to identify cells that derive from the same original clonally-expanded cell. 
For more information on BraCeR, see our [paper in Nature Methods](https://www.nature.com/articles/s41592-018-0082-3) or the [bioRxiv preprint](http://www.biorxiv.org/content/early/2017/09/07/185504) that preceeded it.

BraCeR builds on the well-verified tool for reconstruction of T cell receptor sequences from single cell RNA-seq data (TraCeR).
For more information on TraCeR, its validation and how it can be applied to investigate T cell populations during infection, see our [paper in Nature Methods](http://www.nature.com/nmeth/journal/vaop/ncurrent/full/nmeth.3800.html) or the [bioRxiv preprint](http://biorxiv.org/content/early/2015/08/28/025676) that preceded it.

Please email questions / problems to il5@sanger.ac.uk

## Installation
BraCeR is written in Python and so can just be downloaded, made executable (with `chmod u+x bracer`) and run or run with `python bracer`. Download the latest version and accompanying files from www.github.com/teichlab/bracer. 

BraCeR relies on several additional tools and Python modules that you should install. BraCeR also requires R (>= 3.1.2) and some R packages for lineage reconstruction (optional).

### Pre-requisites

#### OS and hardware requirements
For optimal performance, run BraCeR `assemble` for multiple cells in parallel on a High Performance Computing Cluster with 6+ GB of RAM and 4+ cores. 

The developmental version of BraCeR is tested on Ubuntu 12.04.5 LTS. The Docker version of BraCeR has been tested on Mac OSX, Windows and Linux operating systems.

#### Software requirements
1. [Python3](https://www.python.org) - BraCeR requires Python (>=3.4.0), as one of the required tools has this as a requirement.
2. [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) - required for alignment of reads to synthetic BCR genomes. Bowtie1 is also required.
3. [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - required for assembly of reads into BCR contigs. BraCeR requires Trinity v2.4.0.
4. [IgBLAST](http://www.ncbi.nlm.nih.gov/igblast/faq.html#standalone) - required for analysis of assembled contigs. (ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/).
5. [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) - required for determination of isotype. (ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/).
6. [Kallisto](http://pachterlab.github.io/kallisto/) - software for quantification of BCR expression.
7. [Graphviz](http://www.graphviz.org) - Dot and Neato drawing programs required for visualisation of clonotype graphs. This is optional - see the [`--no_networks` option](#options-1) to [`summarise`](#summarise-summary-and-clonotype-networks).
8. [PHYLIP](http://evolution.genetics.washington.edu/phylip.html) - dnapars program of PHYLIP is required for lineage reconstruction. This is optional - see the [`--infer_lineage` option](#options-1) to [`summarise`](#summarise-summary-and-clonotype-networks).   
9. [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) - required for adapter and quality trimming (optional).

##### Software versions
BraCeR has been tested on the following versions of software dependencies: bowtie2 v2.2.8, bowtie v.1.1.2, IgBlast v.1.4.0 - v.1.7.0, BLAST v.2.2.31+, Kallisto v.0.43.0, Trinity v.2.4.0, graphwiz v.2.26.3, changeo v.0.3.7, RScript v.3.3.2, phylip (dnapars) v.3.696, Trim Galore v.0.4.4, cutadapt v.1.14, ggplot2 v.2.2.1, alakazam v.0.2.6
	

##### Installing IgBlast 
You should also ensure to set the `$IGDATA` environment variable to point to the location of the IgBlast executable. For example run `export IGDATA=/<path_to_igblast>/igblast/1.4.0/bin`.

#### R packages
The following R packages are required if BraCeR is run with `--infer_lineage`.

1. [ggplot2](https://cran.rstudio.com/web/packages/ggplot2/index.html)
2. [Rscript](https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html)
3. [Alakazam](http://alakazam.readthedocs.io/)

## Setup 
To set up the python dependencies, use the requirements file:

    pip3 install -r requirements.txt

It is **highly** recommended that numpy and biopython are first installed through your system's package manager or conda. Setting up the python requirements takes approximately 90 s on a "normal" desktop computer with a "normal" internet connection. 

Note: Seaborn depends on the module statsmodels, which if updated through other packages may cause problems in Seaborn. If such issues arise, try to uninstall statsmodels and install again:
  
    conda uninstall statsmodels --yes
    conda install -c taugspurger statsmodels=0.8.0     

If you plan to run BraCeR with `--infer_lineage` to create lineage trees, please make sure that you have installed R (>= 3.1.2), ggplot2 (>= 2.0.0), Rscript (>=3.3.2) and Alakazam (>= 0.2.7). 

The bracer module is then installed using:

    python setup.py install

This will add the binary 'bracer' to your local bin folder, which can then be run from anywhere. Installing the bracer module takes about a second.

If you would like to contribute to BraCeR, you can set up a development version with

    python setup.py develop

Which will make BraCeR accessible in your python environment, and incorporate local updates to the code.

Once the prerequisites above are installed and working you're ready to tell BraCeR where to find them.

BraCeR uses a configuration file to point it to the locations of files that it needs and a couple of other options.
An example configuration file is included in the repository - `bracer.conf`.
By default, this is `~/.bracerrc`. If bracer fails to find this file, it will use the `bracer.conf` in the repository.
 The `-c` option to the various bracer modules allows you to specify any other file to act as the configuration file.

**Important:** Make sure to edit the configuration file before using BraCeR.

### External tool locations 
BraCeR will look in your system's `PATH` for external tools. You can override this behaviour by editing your `~/.bracerrc`.
Edit `~/.bracerrc` (or a copy) so that the paths within the `[tool_locations]` section point to the executables for all of the required tools.

	[tool_locations]
	#paths to tools used by BraCeR for alignment, quantitation, etc
	bowtie2_path = /path/to/bowtie2
	igblastn_path = /path/to/igblastn
	blastn_path = /path/to/blastn
	kallisto_path = /path/to/kallisto
	trinity_path = /path/to/trinity
	dot_path = /path/to/dot
	neato_path = /path/to/neato
	changeo_path = /path/to/directory_containing_changeo_scripts
	rscript_path = /path/to/Rscript
	dnapars_path = /path/to/dnapars
	trim_galore_path = /path/to/trim_galore
	cutadapt_path = /path/to/cutadapt
		
Make sure that `changeo_path` points to the directory containing the Change-O scripts (`DefineClones.py`, `CreateGermlines.py` and  `MakeDb.py`). 

#### Trinity options 
##### Jellyfish memory 
	[trinity_options]
	#line below specifies maximum memory for Trinity Jellyfish component. Set it appropriately for your environment.
	max_jellyfish_memory = 1G

Trinity needs to know the maximum memory available to it for the Jellyfish component. Specify this here.
  

#### Base transcriptomes for Kallisto 
	[kallisto_transcriptomes]
	Mmus = /path/to/kallisto/transcriptome_for_Mmus
	Hsap = /path/to/kallisto/transcriptome_for_Hsap

Location of the transcriptome fasta file to which the specific BCR sequences will be appended from each cell. This must be a plain-text fasta file so decompress it if necessary. Transcriptome files for human or mice may be downloaded with the following code:

    mkdir GRCh38 && cd GRCh38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz && \
    gunzip gencode.v27.transcripts.fa.gz && python3 /path/to/bracer/docker_helper_files/gencode_parse.py gencode.v27.transcripts.fa && rm gencode.v27.transcripts.fa

    mkdir GRCm38 && cd GRCm38 && wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_mouse/release_M15/gencode.vM15.transcripts.fa.gz && \
    gunzip gencode.vM15.transcripts.fa.gz && python3 /bracer/docker_helper_files/gencode_parse.py gencode.vM15.transcripts.fa && rm gencode.vM15.transcripts.fa

### BraCeR directory

	[bracer_location]
	#Path to where BraCeR was originally downloaded
	bracer_path = /path/to/bracer

Location of the cloned BraCeR repository containing bracerlib, test_data, resources etc. Eg. `/user/software/bracer` if you cloned the bracer repository into `/user/software`.
Needed for localisation of resources and test_data if running BraCeR with the bracer binary.

## Testing BraCeR 
BraCeR comes with a small dataset in `test_data/` (containing only BCR reads for a single cell) that you can use to test your installation and config file and confirm that all the prerequisites are working. Run it as:

    bracer test -p <ncores> -c <config_file>
    
**Note:** The data used in the test are derived from human B cells so make sure that the config file points to the appropriate human base transcriptome file.

You can also pass the following options to change the Graphviz output format, to prevent attempts to draw network graphs or to test BraCeR with lineage reconstruction.

* `-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
* `--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.
* `--infer_lineage` : Run BraCeR with lineage reconstruction.

    
Running `test` will peform the [`assemble`](#assemble-bcr-reconstruction) step using the small test dataset. It will then perform [`summarise`](#summarise-summary-and-clonotype-networks) using the assemblies that are generated along with pre-calculated output for two other cells (in `test_data/results`). Running the test should take approximately 34 minutes on a "normal" desktop computer.

Compare the output in `test_data/results/filtered_BCR_summary` with the expected results in `test_data/expected_summary`. There should be three cells, two with one productive heavy and one productive lambda, and one cell with one productive heavy, one productive lambda and one non-productive kappa. Cells 2 and 3 should be in a clonotype.


## Using BraCeR 
BraCeR has three modes: *assemble*, *summarise* and *build*. 

*Assemble* takes fastq files of paired-end RNA-seq reads from a single-cell and reconstructs BCR sequences.

*Summarise* takes a set of directories containing output from the *assemble* phase (each directory represents a single cell) and summarises BCR recovery rates as well as generating clonotype networks. 

*Build* creates new combinatorial recombinomes for species other than the inbuilt Human and Mouse.


### *Assemble*: BCR reconstruction 

#### Usage 

    bracer assemble [options]  <cell_name> <output_directory> [<file_1>] [<file_2>]


##### Main arguments
* `<cell_name>` : Name of the cell. This is arbitrary text that will be used for all subsequent references to the cell in filenames/labels etc.     
* `<output_directory>` : Directory for output. Will be created if it doesn't exist. Cell-specific output will go into `/<output_directory>/<cell_name>`. This path should be the same for every cell that you want to summarise together.
* `<file_1>` : FASTQ file containing #1 mates from paired-end sequencing or all reads from single-end sequencing.    
* `<file_2>` : FASTQ file containing #2 mates from paired-end sequencing. Do not use if your data are from single-end sequencing.  

##### Options 

* `-p/--ncores <int>` : Number of processor cores available. This is passed to Bowtie2, Trinity, and Kallisto. Default=1.
* `-c/--config_file <conf_file>` : Config file to use. Default = `~/.bracerrc`
* `--resource_dir <resource_dir>`:  The directory containing the resources required for alignment. By default this is the  resources  directory in this repository, but can be pointed to a user-built set of resources.
* `-s/--species` <species> : Species from which the B cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. Default = `Hsap`.
* `-r/--resume_with_existing_files` : If this is set, BraCeR will look for existing output files and not re-run steps that already appear to have been completed. This saves time if BraCeR died partway through a step and you want to resume where it left off. 
* `--assembled_file` <fasta_file> : FASTA file containing already assembled sequences for the cell. Providing this file skips the alignment and assembly steps (default: False). If FASTQ file(s) are provided, BraCeR quantifies the provided assembled sequences.  
* `--single_end` : use this option if your data are single-end reads. If this option is set you must specify fragment length and fragment sd as below.
* `--fragment_length <int>` : Estimated average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--fragment_sd <int>` : Estimated standard deviation of average fragment length in the sequencing library. Used for Kallisto quantification. Required for single-end data. Can also be set for paired-end data if you don't want Kallisto to estimate it directly.
* `--loci`: Space-separated list of loci to reconstruct (default: ['H', 'K', 'L']).
* `--max_junc_len` <int>:  Maximum permitted length of CDR3 nucleotide sequence or junction string. Used to filter out artefacts (default: 100). 
* `--no_trimming` : Skip adapter and quality trimming of raw reads.
* `--keep_trimmed_reads` : Do not delete files containing trimmed raw reads.
 

#### Output 

For each cell, an `/<output_directory>/<cell_name>` directory will be created. This will contain the following subdirectories.

1. `<output_directory>/<cell_name>/trimmed_reads`  
    Subdirectory containing the output from Trim Galore! if assemble is run with `--keep_trimmed_reads`.

2. `<output_directory>/<cell_name>/aligned_reads`  
    This contains the output from Bowtie2 with the sequences of the reads that aligned to the synthetic genomes.

3. `<output_directory>/<cell_name>/Trinity_output`  
    Contains fasta files for each locus where contigs could be assembled. Also two text files that log successful and unsuccessful assemblies.

4. `<output_directory>/<cell_name>/IgBLAST_output`  
    Files with the output from IgBLAST for the contigs from each locus. 

5. `<output_directory>/<cell_name>/BLAST_output`  
    Files with the output from BLAST for the contigs from each locus.

6. `<output_directory>/<cell_name>/unfiltered_BCR_seqs`  
    Files describing the BCR sequences that were assembled prior to filtering by expression if necessary.
    - `unfiltered_BCRs.txt` : text file containing BCR details. Begins with count of productive/total rearrangements detected for each locus. Then details of each detected recombinant.
    - `<cell_name>_BCRseqs.fa` : FASTA file containing reconstructed BCR sequences.
    - `<cell_name>.pkl` : Python [pickle](https://docs.python.org/2/library/pickle.html) file containing the internal representation of the cell and its recombinants as used by BraCeR. This is used in the summarisation steps.

7. `<output_directory>/<cell_name>/expression_quantification`  
    Contains Kallisto output with expression quantification of the entire transcriptome *including* the reconstructed BCRs.

8. `<output_directory>/<cell_name>/filtered_BCR_seqs`  
    Contains the same files as the unfiltered directory above but these recombinants have been filtered so that only the two most highly expressed from each locus are retained. This resolves biologically implausible situtations where more than two recombinants are detected for a locus. **This directory contains the final output with high-confidence BCR assignments**.


### *Summarise*: Summary and clonotype networks 

#### Usage 
    bracer summarise [options] <input_dir>

##### Main argument 
* `<input_dir>` : directory containing subdirectories of each cell you want to summarise. 

##### Options 
* `-c/--config_file <conf_file>` : config file to use. Default = `~/.bracerrc`
* `--resource_dir <resource_dir>` : The directory containing the resources required for alignment. By default this is the  resources  directory in this repository, but can be pointed to a user-built set of resources.
* `-u/--use_unfiltered` : Set this flag to use unfiltered recombinants for summary and networks rather than the recombinants filtered by expression level.   
* `-s/--species` <species> : Species from which the B cells were derived. Options are `Mmus` or `Hsap` for mouse or human data. Default = `Hsap`.
* `--loci`: Space-separated list of loci to summarise (default: ['H', 'K', 'L'])
* `-g/--graph_format` : Output format for the clonotype networks. This is passed directly to Graphviz and so must be one of the options detailed at http://www.graphviz.org/doc/info/output.html.  
* `--no_networks` : Don't try to draw clonotype network graphs. This is useful if you don't have a working installation of Graphviz.
* `--IGH_networks` : Base clonality only on heavy chain (allow different or no light chains in a clonotype). Default=False.
* `--dist <float>` : Distance value to use for clonal inference by Change-O. Heavily mutated datasets may require a higher distance value, whereas datasets enriched in naive B cells may require a lower value. Default=0.2.
* `--include_multiplets` : Do not exclude cells containing more than two recombinants for a locus from downstream analyses, including networks and clonotype analysis. Default=False.
* `--infer_lineage` : Construct lineage trees for clone groups shown in clonal network. Do not use if you do not have dnapars (PHYLIP), R and the required R packages installed. Default=False.


#### Output 
Output is written to `<input_dir>/filtered_BCR_summary` or `<input_dir>/unfiltered_BCR_summary` depending on whether the `--use_unfiltered` option was set.

The following output files and subdirectories may be generated (depending on which arguments BraCeR is run with):

1. `BCR_summary.txt`    
    Summary statistics describing successful BCR reconstruction rates and the numbers of cells with 0, 1, 2 or more recombinants for each locus.
2. `changeodb.tab`    
    Tab-delimited database file containing all reconstructed sequences (except suspected multiplets unless run with `--include_multiplets`) 
3. `filtered_multiplets_changeodb.tab`   
    Tab-delimited database file containing all reconstructed sequences from suspected multiplets (unless run with `--include_multiplets`)
4. `IMGT_gapped.tab`    
    Tab-delimited database file containing information parsed from IgBlast with IMGT-gapped reference sequences for all reconstructed sequences.
5. `reconstructed_lengths_BCR[H|K|L].pdf` and  `reconstructed_lengths_BCR[H|K|L].txt`    
    Distribution plots (and text files with underlying data) showing the lengths of the VDJ regions from assembled BCR contigs. Longer contigs give higher-confidence segment assignments. Text files are only generated if at least one BCR is found for a locus. Plots are only generated if at least two BCRs are found for a locus. 
6. `clonotype_sizes.pdf` and `clonotype_sizes.txt`    
    Distribution of clonotype sizes as bar graph and text file.
7.  `clonotype_network_[with|without]_identifiers.<graph_format>`    
    graphical representation of clonotype networks either with full recombinant identifiers or just lines indicating presence/absence of recombinants.
8.  `clonotype_network_[with|without]_identifiers.dot`    
    files describing the clonotype networks in the [Graphviz DOT language](http://www.graphviz.org/doc/info/lang.html)
9. `lineage_trees/`    
    Subdirectory containing lineage tree output files if run with `--infer_lineage`
10. Intermediate output files from the various steps.
    - `changeo_input_[H|K|L].tab` : Tab-delimited database file used as input for Change-O DefineClones. Contains all productive sequences for the locus.
    - `changeo_input_[H|K|L]_clone-pass.tab` : Output of Change-O DefineClones bygroup.
    - `igblast_input_[H|K|L].fa` : FASTA file containing all sequences that are shared within one of the clone groups in the clonal network. 
    - `igblast_[H|K|L].fmt7` : IgBlast output file for `igblast_input_[H|K|L].fa`, having run IgBlast with IMGT-gapped reference sequences.  
    - `igblast_[H|K|L]_db-modified.tab` : Tab-delimited database file created from `igblast_[H|K|L].fmt7` through Change-O. Modified to include clone group, isotype and cell name columns.  
    - `igblast_[H|K|L]_db-modified_germ-pass.tab` : Output file from germline reconstruction step (running Change-O CreateGermlines on `igblast_[H|K|L]_db-modified.tab`).  
    - `concatenated_lineage_input.tab` : Tab-delimited database file used as input for lineage reconstruction with Alakazam. Contains concatenated heavy and light chain sequences and their inferred germline sequences for each clone group in the clonal networks.

### *Build*: Build Combinatorial Recombinomes for a Given Species

#### Usage
    bracer build <species> <locus_name> <N_padding> <colour> <V_seqs> <J_seqs> <C_seqs> <D_seqs>

#### Main Arguments
* `<species>` : Species (e.g. Hsap).
* `<locus_name>` : Name of locus (e.g. H)
* `<N_padding>` : Number of ambiguous N nucleotides between V and J
* `<colour>` : Colour for the productive recombinants (optional). Specify as HTML (e.g. E41A1C) or use "random"
* `<V_seqs>` : Fasta file containing V gene sequences
* `<J_seqs>` : Fasta file containing J gene sequences
* `<C_seqs>` : Fasta file containing single constant region sequence
* `<D_seqs>` : Fasta file containing D gene sequences (required for heavy chain)

#### Options

* `-f/--force_overwrite` : Force overwrite of existing resources
* `-c/--config_file <conf_file>` : config file to use. Default = `~/.bracerrc`
* `--C_db <alt_C_seqs>` : Specify alternative FASTA file (if other than the one used to make recombinomes) containing all C gene sequences for creation of BLAST database to correctly identify isotype (optional)
* `--V_gapped <gapped_V_seqs>` :  FASTA file containing IMGT-gapped V reference  sequences (optional, but HIGHLY recommended). Required for reliable detection of CDR3, lineage reconstruction and creation of IMGT-gapped tab-delimited databases
* `--igblast_aux <igblast_auxiliary_file>` : IgBlast auxiliary file for species. HIGHLY recommended, required for reliable detection of CDR3, lineage reconstruction and creation of IMGT-gapped tab-delimited databases.
* `--resource_dir <resource_dir>` : Root directory for resources

## Docker image

BraCeR is also available as a standalone Docker image on [DockerHub](https://hub.docker.com/r/teichlab/bracer/), with all of its dependencies installed and configured appropriately. Running BraCeR from the image is very similar to running it from a normal installation. You can pass all the appropriate arguments to the Docker command with the usual syntax as described above. One difference is that you don't need to worry about specifying a configuration file. This is included in the container.

To run the BraCeR Docker image, run the following command from within a directory that contains your input data:

	docker run -it --rm -v $PWD:/scratch -w /scratch teichlab/bracer
	
followed by any arguments that you want to pass to BraCeR. 

The `-it` flag ensures that you see all the information BraCeR prints to the screen during its run, while `--rm` deletes the individual container created based on the image for the purpose of the run once the analysis is complete, not littering your computer's drive. `-v` creates a volume, allowing the created container to see the contents of your current working directory, and the `-w` flag sets the container's working directory to the newly created volume. 

For example, if you wanted to run the test analysis, you should clone this GitHub repository, navigate to its main directory so you can see the `test_data` folder, and call the following (you need to specify the `-o test_data` so that the results get written to the volume you created, ensuring you can see them after the analysis is finished):

	docker run -it --rm -v $PWD:/scratch -w /scratch teichlab/bracer test -o test_data

If you wish to use `bracer build`, you will need to specify `--resource_dir /scratch`, as otherwise the resulting resources will be saved in the default location of the container and subsequently get forgotten about when the build analysis completes, making them unuseable for any actual analyses you may want to perform. This will make the Docker container save the resulting resources in the volume you created, and you can use them for assemble/summarise by running the Dockerised BraCeR from the same directory as the one you used for the build and specifying `--resource_dir /scratch`.

You may need to explicitly tell Docker to increase the memory that it can use. Instructions for [Windows](https://docs.docker.com/docker-for-windows/#advanced) and [Mac](https://docs.docker.com/docker-for-mac/#advanced). Something like 6 or 8 GB is likely to be ok.
