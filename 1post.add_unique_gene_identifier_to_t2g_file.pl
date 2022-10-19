use strict;
use warnings;

# REQUIREMENTS
my $species = 'Rattus norvegicus'; my $species_spaceless = $species; $species_spaceless =~ s/\s/\_/;
my $tx2gene = "indexes/$species_spaceless.transcript_to_gene_lookup.tsv"; # from 1.make_protein_coding_transcriptome.pl
if (!(-e($tx2gene))) { print "ERROR: cannot find $tx2gene\n"; exit 1; }

# OUTPUT
my $out_file = "indexes/$species_spaceless.transcript_to_gene_lookup_with_unique_IDs.tsv";
open(OUT,'>',$out_file) or die $!;
print OUT "Gene name\tGene ID\tTranscript (Ensembl or RefSeq) ID\tDescription\tLocation\tUnique gene identifier\n";

# STORE GENE NAMES FOR EACH TRANSCRIPT ID
my %transcript_to_gene_name_lookup = (); my %gene_name_to_transcript_id_lookup = (); my %gene_name_to_gene_id_lookup = (); my %gene_descs = (); my %gene_coords = ();
# is this gene name assigned to multiple gene IDs? we'll find out by looking at the coordinates to see which ones overlap
my %gene_ids_per_gene_name = ();
open(IN,$tx2gene) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $gene_id = $line[1]; my $loc = $line[4];
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
open(IN,$tx2gene) or die $!;
while(<IN>)
	{ next if ($. == 1);
	  my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $gene_name = $line[0]; my $gene_id = $line[1]; my $loc = $line[4];
	
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
	  my $unique_gene_name = $gene_name;
	  if ($number_of_gene_ids_per_gene_name > 1)
		{ my $new_gene_name = "$gene_name/$gene_id";
		  $unique_gene_name = $new_gene_name;
		}
	  print OUT "$line\t$unique_gene_name\n";
	}
close(IN) or die $!;
close(OUT) or die $!;
exit 1;