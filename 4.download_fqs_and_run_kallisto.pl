=head
AFTER USAGE, WE CAN SUBMIT ALL SCRIPTS TO A SLURM CLUSTER AS SO:

cd sh
for i in $(find sh -mindepth 1 -maxdepth 1 ;); do sbatch $i; done

=cut

use strict;
use warnings;

# REQUIREMENTS
my $homedir  = 'expr_atlas';
my $species  = 'Rattus_norvegicus';
my $metadata = "$homedir/$species.metadata.tsv"; 			  # created by 3.parse_metadata.pl
my $index    = "$homedir/indexes/$species.protein_coding.idx"; # created by 1.make_protein_coding_transcriptome.pl
my $se_len   = "$homedir/determine_read_length_from_fastp.pl"; # another script; used to parse fastp output and report an average read length. This is a necessary workaround for the particular cluster used in this study as 'jq' (which would have allowed the task to be performed in one line) was not available.
my $ssh_file = "$homedir/.aspera/connect/etc/asperaweb_id_dsa.openssh";
my $fatal    = 0;
if (!(-e($metadata))) { $fatal++; print "ERROR: cannot find $metadata\n"; }
if (!(-e($index))) 	  { $fatal++; print "ERROR: cannot find $index\n"; 	  }
if (!(-e($se_len)))   { $fatal++; print "ERROR: cannot find $se_len\n";   }
if (!(-e($ssh_file))) { $fatal++; print "ERROR: cannot find $ssh_file\n"; }
exit 1 if ($fatal > 0);

# PARAMETERS
my $num_seeds     = 5;
my $downsample_to = 10000000;
my $num_threads   = 10;

# OUTPUT
my $sh_dir  = "$homedir/sh_scripts";
my $out_dir = "$homedir/kallisto_output";
if (!(-d($sh_dir)))  			{ mkdir $sh_dir  			or die $!; }
if (!(-d($out_dir))) 			{ mkdir $out_dir 			or die $!; }
if (!(-d("$out_dir/$species"))) { mkdir "$out_dir/$species" or die $!; }

# STORE ALL RUN IDs PER SAMPLE
# if ENA fastq URLs are available, we will download each sample's fqs using aspera
# else if .sra URLs are available, we will download each sample's .sra file using wget and then decompile into fqs using the SRA toolkit
my %samples = (); my %library_layouts = ();
open(IN,$metadata) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $sample_id = $line[1]; my $run_id = $line[2]; my $library_layout = $line[4]; my $fq_url1 = $line[5]; my $fq_url2 = $line[6]; my $sra_url = $line[7];
	  push(@{$samples{$sample_id}},[$run_id,$fq_url1,$fq_url2,$sra_url]);
	  $library_layouts{$sample_id}{$library_layout}++;
	}
close(IN) or die $!;

# FOR EACH SAMPLE...
my @samples = ();
while((my $sample_id,my $irrel)=each(%samples))
	{ push(@samples,$sample_id); }
my @sorted_samples = sort {$a cmp $b} @samples;
for(my $x=0;$x<@sorted_samples;$x++)
	{ next if ($x > 10);
	  print "$x of $#sorted_samples\n";
	  my $sample_id = $sorted_samples[$x];
	  my $num_library_layouts = scalar keys %{$library_layouts{$sample_id}};
	  next if ($num_library_layouts != 1); # CHECKPOINT: skip if the sample ID is associated with both paired-end and single-end fqs; these have already been pre-processed in some way
	  my $library_layout = '';
	  while((my $this_library_layout,my $irrel)=each(%{$library_layouts{$sample_id}}))
		{ $library_layout = $this_library_layout; }
	  
	  # CHECK THAT WE HAVE NOT PROCESSED THIS SAMPLE ALREADY. IF WE HAVE, WE'LL HAVE PRODUCED x KALLISTO OUTPUT FILES, ONE PER SEED.
	  my $already_done = 0;
	  if (-d("$out_dir/$species/$sample_id"))
		{ opendir(DIR,"$out_dir/$species/$sample_id") or die $!;
		  my @files = readdir(DIR);
		  closedir(DIR) or die $!;
		  foreach my $file (@files)
			{ next if (($file eq '.') or ($file eq '..'));
			  if ($file =~ /^$sample_id\.\d+\.tsv$/)
				{ $already_done++; }
			}
		}
	  next if ($already_done == $num_seeds);
	  
	  # CREATE A SLURM-SUITABLE SHELL SCRIPT TO DOWNLOAD FQs AND RUN x ITERATIONS OF KALLISTO
	  open(SH,'>',"$sh_dir/$sample_id.sh") or die $!;
	  print SH "#!/bin/bash\n";
	  print SH "#SBATCH --partition=batch\n";
	  print SH "#SBATCH --job-name=$sample_id\n";
	  print SH "#SBATCH --ntasks=1\n";
	  print SH "#SBATCH --mem=10G\n";
	  print SH "#SBATCH --time=02:00:00\n";
	  print SH "#SBATCH --output=%j_%x.out\n";
	  print SH "#SBATCH --error=%j_%x.err\n\n";
	  print SH "module add aspera/3.9.8\n";
	  print SH "module add seqtk/20201102\n";
	  print SH "module add fastp/0.20.1\n";
	  print SH "module add kallisto/0.46.1\n";
	  print SH "module add sratoolkit/2.9.6\n";
	  print SH "mkdir $out_dir/$species/$sample_id\n";
	  print SH "cd $out_dir/$species/$sample_id\n";
	  
	  # DOWNLOAD ALL FQs FOR EACH RUN OF THIS SAMPLE
	  my $file_line1 = ''; my $file_line2 = '';
	  my @fqs_to_delete = ();
	  my @arr = @{$samples{$sample_id}};
	  for(my $x=0;$x<@arr;$x++)
		{ my $run_id = $arr[$x][0]; my $url1 = $arr[$x][1]; my $url2 = $arr[$x][2]; my $url3 = $arr[$x][3];
		  
		  if ($url1 ne 'NA') # if ENA fastq URLs are available (which they will be if $url1, at minimum, is available), we will download each sample's fqs using aspera
			{ my $fq_1 = $url1; my $fq_2 = $url2;
			  if ($fq_1 =~ /^.+\/(.*?)$/) { $fq_1 = $1; }
			  if ($fq_2 =~ /^.+\/(.*?)$/) { $fq_2 = $1; }
			  print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url1 $fq_1\n";
			  if (($url2 eq 'NA') and ($library_layout eq 'SINGLE'))
				{ $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  push(@fqs_to_delete,$fq_1);
				}
			  elsif (($url2 ne 'NA') and ($library_layout eq 'PAIRED'))
				{ print SH "ascp -QT -l 100m -P33001 -i $ssh_file $url2 $fq_2\n";
				  $file_line1 .= "$out_dir/$species/$sample_id/$fq_1 ";
				  $file_line2 .= "$out_dir/$species/$sample_id/$fq_2 ";
				  push(@fqs_to_delete,$fq_1,$fq_2);
				}
			}
		  elsif ($url3 ne 'NA') # else if .sra URLs are available, we will download each sample's .sra file using wget and then decompile into fqs using the SRA toolkit
			{ print SH "wget $url3 -O $run_id.sra\n";
			  print SH "fastq-dump -I --split-files --gzip $run_id.sra\n";
			  print SH "rm $run_id.sra\n";
			  if ($library_layout eq 'SINGLE')
				{ my $fq = "$run_id"."_1.fastq.gz ";
				  $file_line1 .= "$fq ";
				  push(@fqs_to_delete,$fq);
				}
			  elsif ($library_layout eq 'PAIRED')
				{ my $fq_1 = "$run_id"."_1.fastq.gz ";
				  my $fq_2 = "$run_id"."_2.fastq.gz ";
				  $file_line1 .= "$fq_1 ";
				  $file_line2 .= "$fq_2 ";
				  push(@fqs_to_delete,$fq_1,$fq_2);
				}
			}
		}
	  $file_line1 =~ s/\s+$//; $file_line2 =~ s/\s+$//;
		
	  # IF THERE HAVE BEEN MULTIPLE SETS OF RUN FQS PER SAMPLE, THEN MERGE THEM ALL TOGETHER INTO ONE FINAL SAMPLE.FQ
	  my $fq_1 = "$sample_id.1.fq.gz";
	  my $fq_2 = "$sample_id.2.fq.gz";
	  my $fq_s = "$sample_id.fq.gz";
	  if ($#arr > 0)
		{ if ($library_layout eq 'PAIRED')
			{ print SH "cat $file_line1 > $fq_1\n";
			  print SH "cat $file_line2 > $fq_2\n";
			}
		  elsif ($library_layout eq 'SINGLE')
			{ print SH "cat $file_line1 > $fq_s\n";
			}
		  foreach my $fq (@fqs_to_delete)
			{ print SH "rm $fq\n"; }
		}
	  elsif ($#arr == 0) # ELSE JUST RENAME THE SINGLE RUN.FQ TO SAMPLE.FQ
		{ if ($library_layout eq 'PAIRED')
			{ print SH "mv $file_line1 $fq_1\n";
			  print SH "mv $file_line2 $fq_2\n";
			}
		  elsif ($library_layout eq 'SINGLE')
			{ print SH "mv $file_line1 $fq_s\n";
			}
		}
	  
	  # DOWNSAMPLE THE FQs, x TIMES, AND THEN QUANTIFY EXPRESSION
	  my $max_seeds = $num_seeds-$already_done;
	  if ($library_layout eq 'PAIRED')
		{ for(my $seed_num=1;$seed_num<=$max_seeds;$seed_num++)
			{ my $seed = int(rand(10000000));
			  my $seed_fq1 = "$sample_id.downsampled.1.fq.gz";
			  my $seed_fq2 = "$sample_id.downsampled.2.fq.gz";
			  
			  # downsample fqs
			  print SH "seqtk sample -s $seed $fq_1 10000000 > $seed_fq1\n";
			  print SH "seqtk sample -s $seed $fq_2 10000000 > $seed_fq2\n";
			  
			  # quantify expression
			  print SH "kallisto quant --index=$index --bias --threads $num_threads --output-dir=$sample_id.$seed $seed_fq1 $seed_fq2\n";
			  print SH "mv $sample_id.$seed/abundance.tsv $sample_id.$seed.tsv\n";
			  print SH "rm $sample_id.$seed/abundance.h5 $sample_id.$seed/run_info.json\n";
			  print SH "rm -r $sample_id.$seed\n";
			  print SH "rm $seed_fq1 $seed_fq2\n";
			}
		}
	  elsif ($library_layout eq 'SINGLE')
		{ my $clean_fq   = "$sample_id.cleaned.fq.gz";
		  my $fastp_json = "$sample_id.json";
		  my $fastp_html = "$sample_id.html";
		  
		  # we need to do something a little different with single-end data
		  # for paired-end samples, Kallisto estimates the fragment length from the reads so does not need to be user-specified
		  # for single-end samples, fragment length cannot be empirically derived from read mapping and is assumed to follow a truncated Gaussian distribution with user-specified mean and standard deviation
		  # as such, for the single-end libraries, we considered the mean fragment length to be 1.2 × the median read length and the standard deviation to be 0.1 × the mean fragment length
		  # but what is the median read length?
		  # to obtain this, we use the (v. fast) pre-processing tool fastp with default parameters (this performs some basic read cleaning). We then parse fastp's summary JSON to extract the median read length
		  # we could, if we wished, also make use of the cleaned reads in downstream processing but we choose not to - we found there was little substantive difference to expression estimates when doing so
		  
		  print SH "fastp -i $fq_s -o $clean_fq -j $fastp_json -h $fastp_html --thread $num_threads\n";
		  print SH "rm $fastp_html $clean_fq\n";
		  print SH "median_read_length=\$(perl $se_len $fastp_json)\n";
		  print SH "rm $fastp_json\n";
		  print SH "num1=1.2\n";
		  print SH "fragment_length=\$(echo \"\$num1 * \$median_read_length\" | bc -l) \n"; # bash does not support floating-point multiplication so we need to use bc for this. Discussed at https://stackoverflow.com/questions/11039876/multiplication-on-command-line-terminal
		  print SH "num2=0.1\n";
		  print SH "fragment_sd=\$(echo \"\$num2 * \$fragment_length\" | bc -l)\n";
		  for(my $seed_num=1;$seed_num<=$max_seeds;$seed_num++)
			{ my $seed = int(rand(10000000));
			  my $seed_fq = "$sample_id.downsampled.fq.gz";
			  print SH "seqtk sample -s $seed $fq_s 10000000 > $seed_fq\n";
			  print SH "kallisto quant --index=$index --bias --threads $num_threads --single -l \$fragment_length -s \$fragment_sd --output-dir=$sample_id.$seed $seed_fq\n";
			  print SH "mv $sample_id.$seed/abundance.tsv $sample_id.$seed.tsv\n";
			  print SH "rm $sample_id.$seed/abundance.h5 $sample_id.$seed/run_info.json\n";
			  print SH "rm -r $sample_id.$seed\n";
			  print SH "rm $seed_fq\n";
			}
		}

	  # DELETE THE ORIGINAL FQs; WE NO LONGER HAVE A NEED FOR THEM
	  if ($library_layout eq 'PAIRED')
		{ print SH "rm $fq_1 $fq_2\n"; }
	  elsif ($library_layout eq 'SINGLE')
		{ print SH "rm $fq_s\n"; }
	  
	  close(SH) or die $!;
	}
exit 1;