use strict;
use warnings;
use Acme::Tools qw(median);

# REQUIREMENTS
my $species   = 'Rattus_norvegicus';
my $in_dir    = 'kallisto_output'; 								  # created by 4.download_fqs_and_run_kallisto.pl
my $locs_file = "indexes/$species.transcript_to_gene_lookup.tsv"; # created by 1.make_protein_coding_transcriptome.pl
my $fatal     = 0;
if (!(-d($in_dir)))    { $fatal++; print "ERROR: can't find $in_dir\n";    }
if (!(-e($locs_file))) { $fatal++; print "ERROR: can't find $locs_file\n"; }
exit 1 if ($fatal > 0);

# OUTPUT
my $out_file = "$species.avg_expression_per_gene_per_sample.tsv";
open(OUT,'>',$out_file) or die $!;

# STORE GENE NAMES FOR EACH TRANSCRIPT ID
my %transcript_to_gene_name_lookup = (); my %gene_name_to_transcript_id_lookup = (); my %gene_name_to_gene_id_lookup = (); my %gene_descs = (); my %gene_coords = ();
# is this gene name assigned to multiple gene IDs? we'll find out by looking at the coordinates to see which ones overlap
my %gene_ids_per_gene_name = ();
open(IN,$locs_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $gene_id = $line[1]; my $transcript_id = $line[2]; my $desc = $line[3]; my $loc = $line[4];
	  $gene_ids_per_gene_name{$gene_name}{$gene_id} = $loc;
	}
close(IN) or die $!;
my %overlapping_gene_ids = ();
while((my $gene_name,my $irrel)=each(%gene_ids_per_gene_name))
	{ my @gene_locs = ();
	  while((my $gene_id,my $irrel)=each(%{$gene_ids_per_gene_name{$gene_name}}))
		{ my $loc = $gene_ids_per_gene_name{$gene_name}{$gene_id};
		  if ($loc =~ /^(.*?)\:(\d+)\-(\d+):(.*?)$/)
			{ my $chr = $1; my $start = $2; my $end = $3; my $strand = $4;
			  push(@gene_locs,[$gene_id,$chr,$start,$end,$strand]);
			}
		}
	  my @sorted_gene_locs = map { $_->[0] } sort { $a->[1] cmp $b->[1] } map { [$_, $_->[0]] } @gene_locs;
	  for(my $x=0;$x<@sorted_gene_locs;$x++)
		{ my $gene_id = $gene_locs[$x][0]; my $chr = $gene_locs[$x][1]; my $start = $gene_locs[$x][2]; my $end = $gene_locs[$x][3]; my $strand = $gene_locs[$x][4];
		  for(my $y=0;$y<@sorted_gene_locs;$y++)
			{ next if ($x == $y);
			  my $gene_id2 = $gene_locs[$y][0]; my $chr2 = $gene_locs[$y][1]; my $start2 = $gene_locs[$y][2]; my $end2 = $gene_locs[$y][3]; my $strand2 = $gene_locs[$y][4];
			  if (($chr =~ /^$chr2$/) and ($strand == $strand2))
				{ if ( (($start >= $start2) and ($start <= $end2)) or (($end >= $start2) and ($end <= $end2)) or (($start <= $start2) and ($end >= $end2)) or (($start2 >= $start) and ($end2 <= $end)) ) # right-offset, left-offset, X-in-Y, Y-in-X
					{ $overlapping_gene_ids{$gene_name}{$gene_id}{$gene_id2}++;
					  $overlapping_gene_ids{$gene_name}{$gene_id2}{$gene_id}++;
					  #print "$gene_name --> $gene_id ($chr:$start-$end:$strand) overlaps $gene_id2 ($chr2:$start2-$end2:$strand2)\n";
					}
				}
			}
		}
	}
open(IN,$locs_file) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $gene_id = $line[1]; my $transcript_id = $line[2]; my $desc = $line[3]; my $loc = $line[4];
	
	  # does this $gene_id overlap with ONE other, Ensembl, gene ID? it likely will do if $gene_id is NCBI and $gene_id2 is Ensembl
	  my %revised_gene_ids = (); my %discarded_gene_ids = ();
	  while((my $gene_id,my $irrel)=each(%{$overlapping_gene_ids{$gene_name}}))
		{ while((my $gene_id2,my $irrel)=each(%{$overlapping_gene_ids{$gene_name}{$gene_id}}))
			{ if ($gene_id2 =~ /^ENS.*?$/)
				{ $revised_gene_ids{$gene_id}{$gene_id2}++; # i.e. $gene_id is now $gene_id2
				  $discarded_gene_ids{$gene_id}++;
				}
			}
		}
	  my $number_of_revised_gene_ids = scalar keys %{$revised_gene_ids{$gene_id}};
	  my $new_gene_id = $gene_id;
	  if ($number_of_revised_gene_ids == 1)
		{ while((my $gene_id2,my $irrel)=each(%{$revised_gene_ids{$gene_id}}))
			{ $new_gene_id = $gene_id2; }
		}
	  $gene_id = $new_gene_id;
	  
	  my $number_of_gene_ids_per_gene_name = 0;
	  while((my $gene_id,my $irrel)=each(%{$gene_ids_per_gene_name{$gene_name}}))
		{ next if (exists($discarded_gene_ids{$gene_id}));
		  $number_of_gene_ids_per_gene_name++;
		}
	  if ($number_of_gene_ids_per_gene_name > 1)
		{ my $new_gene_name = "$gene_name/$gene_id";
		  $gene_name = $new_gene_name;
		}
		
	  $gene_name_to_transcript_id_lookup{$gene_name}{$transcript_id}++;
	  $transcript_to_gene_name_lookup{$transcript_id} = $gene_name;
	  $gene_name_to_gene_id_lookup{$gene_name} = $gene_id;
	  $gene_descs{$gene_name} = $desc;
	  $gene_coords{$gene_id} = $loc;
	}
close(IN) or die $!;

# STORE EXPRESSION LEVEL DATA
my %tpm_values = (); my %sample_ids = ();
opendir(DIR,"$in_dir/$species") or die $!;
my @sample_ids = readdir(DIR);
closedir(DIR) or die $!;
my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
my $sample_ids_seen = 0; my $sample_ids_total = @sorted_sample_ids; $sample_ids_total = $sample_ids_total-2;
foreach my $sample_id (@sorted_sample_ids)
	{ next if (($sample_id eq '.') or ($sample_id eq '..'));
	  $sample_ids_seen++;
	  next if (!(-d("$in_dir/$species/$sample_id")));
	  print "$species - $sample_id - $sample_ids_seen of $sample_ids_total\n";
	  opendir(DIR,"$in_dir/$species/$sample_id") or die $!;
	  my @files = readdir(DIR);
	  closedir(DIR) or die $!;
      foreach my $sample_id_and_seed (@files)
		{ next if (($sample_id_and_seed eq '.') or ($sample_id_and_seed eq '..'));
		  if ($sample_id_and_seed =~ /^.*?\.(\d+)\.tsv$/)
			{ my $seed = $1;
			  my $abundances_file = "$in_dir/$species/$sample_id/$sample_id.$seed.tsv";
			  next if (!(-e($abundances_file)));
			  open(IN,$abundances_file) or die $!;
			  while(<IN>)
				{ next if ($. == 1);
				  my $line = $_; chomp($line);
				  my @line = split(/\t/,$line);
				  my $transcript_id = $line[0]; my $tpm = $line[4];
				  my $gene_name = '';
				  if (exists($transcript_to_gene_name_lookup{$transcript_id}))
					{ $gene_name = $transcript_to_gene_name_lookup{$transcript_id}; }
				  next if ($gene_name eq '');
				  $tpm_values{$gene_name}{$sample_id}{$seed} += $tpm;
				  $sample_ids{$sample_id}{$seed}++;
				}
			  close(IN) or die $!;
			}
		}
	}
  
# OUTPUT EXPRESSION LEVEL DATA, 1: PRINT HEADERS
my $header_line = "Gene name (if multiple genes have the same name, shown as 'gene name/gene ID')\tGene ID (Ensembl if available, else NCBI)\tTranscript IDs (Ensembl and/or NCBI)\tNo. of transcripts\tNo. of Ensembl-sourced transcripts\tNo. of NCBI-sourced transcripts\tDescription\tCoordinates";
my $out_line = '';
@sample_ids = ();
while((my $sample_id,my $irrel)=each(%sample_ids))
	{ push(@sample_ids,$sample_id); }
@sorted_sample_ids = sort {$a cmp $b} @sample_ids;
foreach my $sample_id (@sorted_sample_ids)
	{ $out_line .= "$sample_id\t"; }
$out_line =~ s/\t$//;
print OUT "$header_line\t$out_line\n";

# OUTPUT EXPRESSION LEVEL DATA, 2: PRINT ALL TPM ESTIMATES
my @gene_names = ();
while((my $gene_name,my $irrel)=each(%tpm_values))
	{ push(@gene_names,$gene_name); }
my @sorted_gene_names = sort {"\L$a" cmp "\L$b"} @gene_names;
foreach my $gene_name (@sorted_gene_names)
	{ my $out_line = '';
	  my @sample_ids = ();
	  while((my $sample_id,my $irrel)=each(%{$tpm_values{$gene_name}}))
		{ push(@sample_ids,$sample_id); }
	  my @sorted_sample_ids = sort {$a cmp $b} @sample_ids;
	  foreach my $sample_id (@sorted_sample_ids)
		{ my @seeds = ();
		  while((my $seed,my $irrel)=each(%{$sample_ids{$sample_id}}))
			{ push(@seeds,$seed); }
		  my @sorted_seeds = sort {$a <=> $b} @seeds;
		  my @tpms = ();
		  foreach my $seed (@sorted_seeds)
			{ my $tpm = 0;
			  if (exists($tpm_values{$gene_name}{$sample_id}{$seed}))
				{ $tpm = $tpm_values{$gene_name}{$sample_id}{$seed}; }
			  if ($tpm !~ /\d+/)
				{ $tpm = 0; }
			  push(@tpms,$tpm);
			}
		  my $avg_tpm = median(@tpms);
		  $avg_tpm = sprintf("%.3f",$avg_tpm);
		  $out_line .= "$avg_tpm\t";
		}
	  $out_line =~ s/\t$//;
	  my $gene_id = $gene_name_to_gene_id_lookup{$gene_name};
	  my $desc    = $gene_descs{$gene_name};
	  my $loc     = $gene_coords{$gene_id};
	  my $no_of_ensembl_transcript_ids = 0; my $no_of_ncbi_transcript_ids = 0;
	  my @transcript_ids = ();
	  while((my $transcript_id,my $irrel)=each(%{$gene_name_to_transcript_id_lookup{$gene_name}}))
		{ push(@transcript_ids,$transcript_id);
		  if ($transcript_id =~ /^ENS.*?$/)
			{ $no_of_ensembl_transcript_ids++; }
		  else
			{ $no_of_ncbi_transcript_ids++; }
		}
	  my @sorted_transcript_ids = sort {$a cmp $b} @transcript_ids;
	  my $no_of_transcript_ids  = @sorted_transcript_ids;
	  my $transcript_ids = join(", ",@sorted_transcript_ids);
	  print OUT "$gene_name\t$gene_id\t$transcript_ids\t$no_of_transcript_ids\t$no_of_ensembl_transcript_ids\t$no_of_ncbi_transcript_ids\t$desc\t$loc\t$out_line\n";
	}
close(OUT) or die $!;
exit 1;