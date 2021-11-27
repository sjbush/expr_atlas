# expr_atlas
**Species non-specific expression atlas creation**

This repository contains the scripts used to create the expression atlas detailed in ["Generation and network analysis of an RNA-seq transcriptional atlas for the rat"](https://www.biorxiv.org/content/10.1101/2021.11.07.467633v1). These scripts are extensible to other species, with older versions previously used to generate expression atlases for [sheep](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1006997), [buffalo](https://www.frontiersin.org/articles/10.3389/fgene.2019.00668/full), [chicken](https://link.springer.com/article/10.1186/s12864-018-4972-7), [goat](https://www.frontiersin.org/articles/10.3389/fgene.2019.01080/full), [pig](https://www.frontiersin.org/articles/10.3389/fgene.2019.01355/full), and the [mononuclear phagocyte system of the mouse](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3000859).
Scripts should be considered "lab quality" rather than "developer quality" and while functional, have not been optimised for speed or memory use. They perform the following steps and have the following software prerequisites, which should be accessible on the command line. More specific detail on prerequisites is available as comments within the scripts themselves.

**Script 1: make_protein_coding_transcriptome.pl**

This takes as input one set of Ensembl cDNAs and one set of NCBI mRNA RefSeqs, from the same annotation (in this example, Rnor_6.0) plus the associated GFFs.
Integrating the two, it produces two output files: (a) a non-redundant set of protein-coding transcripts, to be indexed with [Kallisto](https://pachterlab.github.io/kallisto/about), and (b) a transcript-to-gene lookup table, detailing the sequences within.

**Script 2: parse_bioproject_summary_file.pl**

This script parses the daily-updated [NCBI BioProject ‘summary’ file]( ftp://ftp.ncbi.nlm.nih.gov/bioproject/summary.txt) and, for those BioProjects whose data type is ‘transcriptome or gene expression’ and whose taxonomy ID is, e.g., 10116 (rat), uses [pysradb](https://github.com/saketkc/pysradb) to extract the associated metadata. The output is a directory of metadata files, where available, one per BioProject.

**Script 3: parse_metadata.pl**

This script parses the directory of metadata files to create a summary table of usable samples. This lists, per SRA sample and run ID, the associated URLs for each fq and key metadata categories: tissue/cell type, age, and sex. Samples are retained only if ‘instrument model’ is ‘ILLUMINA’, ‘library source’ is ‘TRANSCRIPTOMIC’, ‘library strategy’ is ‘RNA-Seq’, and the metadata specifies, at minimum, ‘tissue’ or ‘source’. This script can rapidly create a shortlist of RNA-seq libraries for inclusion into an atlas but it is important to note that as there is no standardised nomenclature for NCBI metadata, it is also necessary to obtain metadata by manual review (discussed further in [this article](https://www.frontiersin.org/articles/10.3389/fgene.2019.01355/full)). In this respect, the summary table output by this script should be considered a template for manual refinement.

**Script 4: download_fqs_and_run_kallisto.pl**

This script takes as input the sample metadata table produced by the previous script. It produces a directory of shell scripts, one per sample, which in this case are each formatted for submission to a SLURM cluster. Each script in this batch executes a series of commands per sample: downloading sequencing data (using [Aspera](https://www.ibm.com/aspera/connect/) and the [SRA toolkit](https://github.com/ncbi/sra-tools), where appropriate), pre-processing (using [fastp](https://github.com/OpenGene/fastp)), downsampling (using [seqtk](https://github.com/lh3/seqtk)) to _x_ million reads _y_ times, and finally quantifying expression using [Kallisto](https://pachterlab.github.io/kallisto/about). The output is a directory containing one subdirectory per sample ID, within which are the Kallisto output files (one per downsampled replicate).

**Script 5: parse_kallisto_output.pl**

This script takes as input the directory of Kallisto output files. It produces as output one large table: expression per gene per sample, averaged across all replicates.
