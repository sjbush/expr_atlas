use strict;
use warnings;

# REQUIREMENTS
my $species = 'Rattus_norvegicus';
my $in_dir  = "metadata/$species"; # created by 2.parse_bioproject_summary_file.pl
if (!(-d($in_dir))) { print "ERROR: unable to find $in_dir\n"; exit 1; }

# OUTPUT
my $out_file = "$species.metadata.tsv";
open(OUT,'>',$out_file) or die $!;
print OUT "Sample number\tSRA sample ID\tSRA run ID\tBioProject ID\tLibrary layout\tENA fq URL 1\tENA fq URL 2\tSRA URL\tSource name\tTissue/cell type\tAge\tSex\n";

opendir(DIR,$in_dir) or die $!;
my @files = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_files = sort {$a cmp $b} @files;
my $files_seen = 0; my $files_total = @sorted_files; $files_total = $files_total-2;
my %bioprojects_seen = (); my %samples_seen = ();
foreach my $file (@sorted_files)
	{ next if (($file eq '.') or ($file eq '..'));
	  $files_seen++;
	  print "$files_seen of $files_total\n";
	  my $bioproject_id = '';
	  if ($file =~ /^(.*?)\.txt$/) { $bioproject_id = $1; }
	  my $sample_id_num; my $run_num; my $spot_num; my $library_strategy_num; my $library_source_num; my $library_layout_num; my $platform_num; my $taxid_num; my $download_path0_num; my $download_path1_num; my $download_path2_num; my $download_path3_num;
	  my $source_num; my $tissue_num; my $age_num; my $gender_num; my $sex_num;
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  if ($. == 1)
			{ for(my $x=0;$x<@line;$x++)
				{ my $val = $line[$x];
				  if 	($val eq 'sample_accession')  	  { $sample_id_num		   = $x; }
				  elsif ($val eq 'run_accession')     	  { $run_num 	           = $x; }
				  elsif ($val eq 'total_spots')       	  { $spot_num	           = $x; } # NOTE: spots are not the same as reads, but they're often interchangeable. See https://www.biostars.org/p/12047/ and the SRA glossary (https://www.ncbi.nlm.nih.gov/books/NBK54984/): "reads that are mate pairs are concatenated into a single monolithic 'spot' sequence."
				  elsif ($val eq 'library_strategy')  	  { $library_strategy_num  = $x; }
				  elsif ($val eq 'library_source')	  	  { $library_source_num    = $x; }
				  elsif ($val eq 'library_layout')    	  { $library_layout_num    = $x; }				  
				  elsif ($val eq 'instrument_model_desc') { $platform_num 		   = $x; }
				  elsif ($val eq 'organism_taxid ')	  	  { $taxid_num   		   = $x; } # NOTE: the space in 'organism_taxid ' is NOT an error but is necessary to parse these files
				  elsif ($val eq 'ena_fastq_ftp')         { $download_path0_num    = $x; }
				  elsif ($val eq 'ena_fastq_ftp_1')       { $download_path1_num    = $x; }
				  elsif ($val eq 'ena_fastq_ftp_2')       { $download_path2_num    = $x; }
				  elsif ($val eq 'sra_url')         	  { $download_path3_num    = $x; }
				  elsif ($val eq 'source_name')			  { $source_num			   = $x; }
				  elsif ($val eq 'tissue')			  	  { $tissue_num			   = $x; }
				  elsif ($val eq 'gender')			  	  { $gender_num			   = $x; }
				  elsif ($val eq 'sex')			  	  	  { $sex_num			   = $x; }
				  elsif ($val eq 'age')			  	  	  { $age_num			   = $x; }
				}
			}
		}
	  close(IN) or die $!;
	  if ( (!(defined($sample_id_num))) or (!(defined($run_num))) or (!(defined($spot_num))) or (!(defined($library_strategy_num))) or (!(defined($library_source_num))) or (!(defined($library_layout_num))) or (!(defined($platform_num))) or (!(defined($taxid_num))) or (!(defined($download_path0_num))) or (!(defined($download_path1_num))) or (!(defined($download_path2_num))) or (!(defined($download_path3_num))) )
		{ print "unable to proceed with $bioproject_id: insufficient metadata\n"; }
	  if ( (!(defined($tissue_num))) and (!(defined($source_num))) )
		{ print "unable to process with $bioproject_id: neither 'tissue' nor 'source' stated\n"; }
	  next if ( (!(defined($tissue_num))) and (!(defined($source_num))) ); # CHECKPOINT: discard samples with no recognised tissue/cell origin
	  open(IN,"$in_dir/$file") or die $!;
	  while(<IN>)
		{ next if ($. == 1);
		  my $line = $_; chomp($line);
		  my @line = split(/\t/,$line);
		  next if ( (!(defined($sample_id_num))) or (!(defined($run_num))) or (!(defined($spot_num))) or (!(defined($library_strategy_num))) or (!(defined($library_source_num))) or (!(defined($library_layout_num))) or (!(defined($platform_num))) or (!(defined($taxid_num))) or (!(defined($download_path0_num))) or (!(defined($download_path1_num))) or (!(defined($download_path2_num))) );
		  my $sample_id = $line[$sample_id_num]; my $run_id = $line[$run_num]; my $number_of_reads = $line[$spot_num]; my $library_strategy = $line[$library_strategy_num]; my $library_source = $line[$library_source_num]; my $library_layout = $line[$library_layout_num];
		  my $platform = $line[$platform_num]; my $taxid = $line[$taxid_num];
		  next if ($taxid != 10116); # CHECKPOINT: exclude data not from the rat
		  my $url0 = ''; my $url1 = ''; my $url2 = ''; my $url3 = '';
		  if (defined($line[$download_path0_num])) { $url0 = $line[$download_path0_num]; }
		  if (defined($line[$download_path1_num])) { $url1 = $line[$download_path1_num]; }
		  if (defined($line[$download_path2_num])) { $url2 = $line[$download_path2_num]; }
		  if (defined($line[$download_path3_num])) { $url3 = $line[$download_path3_num]; }
		  my $source = 'unknown/not stated'; my $tissue = 'unknown/not stated'; my $age = 'unknown/not stated'; my $sex = ''; my $gender = '';
		  if   (defined($source_num)) { $source = $line[$source_num];   } # $line[$source_num] may be 'undef', hence the if (!(defined)) below
		  if   (defined($tissue_num)) { $tissue = $line[$tissue_num];   }
		  if   (defined($gender_num)) { $gender = $line[$gender_num];   }
		  if   (defined($sex_num)) 	  { $sex 	= $line[$sex_num];      }
		  if   (defined($age_num)) 	  { $age 	= $line[$age_num];      }
		  if (!(defined($source)))    { $source = 'unknown/not stated'; }
		  if (!(defined($tissue)))    { $tissue = 'unknown/not stated'; }
		  if (!(defined($gender)))    { $gender = 'unknown/not stated'; }
		  if (!(defined($sex)))    	  { $sex    = 'unknown/not stated'; }
		  if (!(defined($age)))       { $age    = 'unknown/not stated'; }
		  if (($sex eq '') and ($gender ne '') and ($gender ne 'unknown/not stated')) { $sex = $gender; }
		  if  ($sex eq '') { $sex = 'unknown/not stated'; }
		  next if (($source eq 'unknown/not stated') and ($tissue eq 'unknown/not stated'));
		  next if ($number_of_reads !~ /\d+/);
		  next if ($platform 		 ne 'ILLUMINA');
		  next if ($library_source   ne 'TRANSCRIPTOMIC');
		  next if ($library_strategy ne 'RNA-Seq');
		  if ($library_layout eq 'SINGLE')
			{ if (($url0 ne '') and ($url1 eq ''))
				{ $url1 = $url0; }
			  $url2 = 'NA';
			}
		  $url1 =~ s/^http:\/\///; $url2 =~ s/^http:\/\///; $url3 =~ s/^http:\/\///;
		  if (($url1 eq '') and ($url2 eq '') and ($url3 eq ''))
			{ print "unable to proceed with $bioproject_id, sample $sample_id, because there are no URLs\n"; }
		  next if (($url1 eq '') and ($url2 eq '') and ($url3 eq '')); # CHECKPOINT: this can only be the case if there's been a layout-labelling error
		  if ($url1 eq '') { $url1 = 'NA'; }
		  if ($url2 eq '') { $url2 = 'NA'; }
		  if ($url3 eq '') { $url3 = 'NA'; }
		  $bioprojects_seen{$bioproject_id}++; $samples_seen{$sample_id}++;
		  my $sample_num = scalar keys %samples_seen;
		  print OUT "$sample_num\t$sample_id\t$run_id\t$bioproject_id\t$library_layout\t$url1\t$url2\t$url3\t$source\t$tissue\t$age\t$sex\n";
		}
	  close(IN) or die $!;
	}
close(OUT) or die $!;

my $num_bioprojects_seen = scalar keys %bioprojects_seen;
my $num_samples_seen 	 = scalar keys %samples_seen;

print "Total BioProjects seen: $num_bioprojects_seen. Total samples seen: $num_samples_seen\n";
my @bioproject_ids = ();
while((my $bioproject_id,my $irrel)=each(%bioprojects_seen))
	{ push(@bioproject_ids,$bioproject_id); }
my @sorted_bioproject_ids = sort {$a cmp $b} @bioproject_ids;
print "BioProjects:\n";
foreach my $bioproject_id (@sorted_bioproject_ids)
	{ print "$bioproject_id\n";
	}

exit 1;