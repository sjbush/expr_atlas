use strict;
use warnings;

# REQUIREMENTS
my $in_file = 'summary.txt'; # FROM: ftp://ftp.ncbi.nlm.nih.gov/bioproject/summary.txt. You get to this file by going to "Download (FTP)" on https://www.ncbi.nlm.nih.gov/bioproject/
if (!(-e($in_file))) { print "ERROR: cannot find $in_file\n"; exit 1; }

# PARAMETERS
my $species = 'Rattus norvegicus';
my $species_spaceless = $species; $species_spaceless =~ s/\s/\_/g;

# OUTPUT
my $out_dir = 'metadata';
if (!(-d($out_dir))) 					  { mkdir  $out_dir 					or die $!; }
if (!(-d("$out_dir/$species_spaceless"))) { mkdir "$out_dir/$species_spaceless" or die $!; }
my $sh_file = 'identify_metadata_for_each_bioproject.sh';
open(SH,'>',$sh_file) or die $!;
print SH "#!/bin/bash\n";

# PARSE THE NCBI SRA BioProject SUMMARY FILE
my $eof;
open(IN,$in_file) or die $!;
while(<IN>) { $eof = $.; }
close(IN) or die $!;
open(IN,$in_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $pc = sprintf("%.4f",(($./$eof)*100)); print "$pc%\n";
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $org = $line[0]; my $ncbi_taxonomy_id = $line[1]; my $bioproject_accession = $line[2]; my $data_type = $line[5];
	  next if (-e("$out_dir/$species_spaceless/$bioproject_accession.txt")); # CHECKPOINT: skip because we've seen this before
	  next if ($data_type ne 'Transcriptome or Gene expression'); # CHECKPOINT: restrict data to only transcriptomic studies (typically this data type refers to RNA-seq, but not always: q-PCR and microarrays may also feature)
	  next if ($ncbi_taxonomy_id != 10116); # CHECKPOINT: restrict data to a species of interest. We do so by reference to NCBI taxonomy ID rather than a simple test of whether $org eq $species because (a) text is more prone to typographical error, and (b) data may be attributed to the names of particular subspecies, e.g. "Canis lupus" vs. "Canis lupus familiaris"
	  print SH "pysradb metadata --detailed $bioproject_accession --saveto $out_dir/$species_spaceless/$bioproject_accession.txt\n";
	}
close(IN) or die $!;
close(SH) or die $!;
exit 1;