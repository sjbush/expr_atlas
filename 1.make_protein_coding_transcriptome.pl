=head
ABOUT THIS SCRIPT:
This script requires one set of Ensembl transcripts (cDNA) and one set of NCBI transcripts (mRNA), from the same annotation (in this case, Rnor_6.0) plus the associated gffs.
Integrating the two, it will produce two output files: (a) a non-redundant set of protein-coding transcripts, and (b) a transcript-to-gene lookup table, detailing the sequences within.

After running this script, a Kallisto index can be built from the combined set of transcripts:
kallisto index -i Rattus_norvegicus.protein_coding.idx Rattus_norvegicus.protein_coding.fa
=cut

use strict;
use warnings;

# REQUIREMENTS
my $species  = 'Rattus norvegicus'; my $species_spaceless = $species; $species_spaceless =~ s/\s/\_/; # we require $species to contain whitespace when parsing $in_file2 (see line 149) but we'll remove the spaces when using this variable in file names
my $in_file1 = "indexes/Rattus_norvegicus.Rnor_6.0.cdna.all.fa"; # FROM: http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz
my $ens_loc  = "indexes/Rattus_norvegicus.Rnor_6.0.104.gff3";    # FROM: http://ftp.ensembl.org/pub/release-104/gff3/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gff3.gz
my $in_file2 = "indexes/GCF_000001895.5_Rnor_6.0_rna.fna"; 		 # FROM: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/all_assembly_versions/suppressed/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_rna.fna.gz (archival versions of the rat genome can be found via https://www.ncbi.nlm.nih.gov/genome?term=rattus%20norvegicus)
my $ncbi_loc = "indexes/GCF_000001895.5_Rnor_6.0_genomic.gff";   # FROM: https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Rattus_norvegicus/all_assembly_versions/suppressed/GCF_000001895.5_Rnor_6.0/GCF_000001895.5_Rnor_6.0_genomic.gff.gz
my $code 	 = 'genetic_code.txt'; 								 # manually created. This file is required when parsing NCBI RefSeqs as we need to remove UTRs from their mRNAs. To do so, we need to first identify the longest ORF.
my $ncbi_to_ens_chr_lookup = 'ncbi_to_ens_chr_lookup.txt'; 		 # manually created. This file is required when integrating data from Ensembl and NCBI as the latter uses a different nomenclature to the former. The table in this file reconciles the two.
if (!(-e($code)))	  { print "ERROR: cannot find $code\n"; 	exit 1; }
if (!(-e($ens_loc)))  { print "ERROR: cannot find $ens_loc\n";  exit 1; }
if (!(-e($ncbi_loc))) { print "ERROR: cannot find $ncbi_loc\n"; exit 1; }
if (!(-e($in_file1))) { print "ERROR: cannot find $in_file1\n"; exit 1; }
if (!(-e($in_file2))) { print "ERROR: cannot find $in_file2\n"; exit 1; }
if (!(-e($ncbi_to_ens_chr_lookup))) { print "ERROR: cannot find $ncbi_to_ens_chr_lookup\n"; exit 1; }

# OUTPUT
my $out_file = "indexes/$species_spaceless.protein_coding.fa";
my $tx2gene  = "indexes/$species_spaceless.transcript_to_gene_lookup.tsv";
open(OUT,'>',$out_file) or die $!; open(TX2GENE,'>',$tx2gene) or die $!;
print TX2GENE "Gene name\tGene ID\tTranscript (Ensembl or RefSeq) ID\tDescription\tLocation\n";

# STORE GENETIC CODE
my %trans = ();
open(CODE,$code) or die "Unable to open $code\n";
while ($_= <CODE>)
	{ /^\s*([ACGTN]{3})\s+([A-Z\*\1-8])/ or die "Error in $code in line $_\n";
      $trans{$1} = $2;
    }
close(CODE) or die $!;

# STORE ENSEMBL GENE COORDINATES
my %ens_gene_locs = ();
open(IN,$ens_loc) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line =~ /^\#/);
	  my @line = split(/\t/,$line);
	  my $chr = $line[0]; my $type = $line[2]; my $start = $line[3]; my $end = $line[4]; my $strand = $line[6]; my $info = $line[8];
	  next if ($type ne 'gene');
	  my $revised_strand = $strand;
	  if 	($strand eq '+') { $revised_strand = 1;  }
	  elsif ($strand eq '-') { $revised_strand = -1; }
	  $strand = $revised_strand;
	  my $gene_id   = ''; if ($info =~ /^.*?ID\=gene\:(.*?)\;.*?$/) { $gene_id   = $1; }
	  my $gene_name = ''; if ($info =~ /^.*?Name\=(.*?)\;.*?$/) 	{ $gene_name = $1; }
	  my $loc = "$chr:$start-$end:$strand";
	  $ens_gene_locs{$gene_id} = $loc;
	  $ens_gene_locs{$gene_name} = $loc;
	}
close(IN) or die $!;

# STORE A LOOKUP TABLE FOR THE NCBI CHROMOSOMES, SO THAT THEY CAN BE SYNCED UP WITH THE ENSEMBL CHROMOSOMES
my %ncbi_to_ens_chr_lookup = ();
open(IN,$ncbi_to_ens_chr_lookup) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  my @line = split(/\t/,$line);
	  my $this_species = $line[0]; my $ens_chr = $line[1]; my $ncbi_chr = $line[2];
	  next if ($this_species ne $species);
	  $ncbi_to_ens_chr_lookup{$ncbi_chr} = $ens_chr;
	}
close(IN) or die $!;

# STORE NCBI GENE COORDINATES AND GENE IDS
my %ncbi_gene_locs = (); my %ncbi_transcript_id_to_gene_id_lookup = ();
open(IN,$ncbi_loc) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  next if ($line =~ /^\#/);
	  my @line = split(/\t/,$line);
	  my $ncbi_chr = $line[0]; my $type = $line[2]; my $start = $line[3]; my $end = $line[4]; my $strand = $line[6]; my $info = $line[8];
	  my $revised_strand = $strand;
	  if 	($strand eq '+') { $revised_strand = 1;  }
	  elsif ($strand eq '-') { $revised_strand = -1; }
	  $strand = $revised_strand;
	  next if (!(exists($ncbi_to_ens_chr_lookup{$ncbi_chr}))); # CHECKPOINT: we are only going to include NCBI genes that map to the Ensembl assembly
	  my $ens_chr = $ncbi_to_ens_chr_lookup{$ncbi_chr};
	  my $loc = "$ens_chr:$start-$end:$strand";
	  if ($type eq 'mRNA')
		{ if ($info =~ /^.*?GeneID\:(.*?)\,Genbank\:(.*?)\..*?$/)
			{ my $gene_id = $1; my $transcript_id = $2;
			  $ncbi_transcript_id_to_gene_id_lookup{$transcript_id} = $gene_id;
			  $ncbi_gene_locs{$gene_id} = $loc;
			}
		}
	}
close(IN) or die $!;

# STORE PROTEIN-CODING TRANSCRIPTS (ENSEMBL)
my %ensembl_seqs = (); my %ensembl_seqs_per_gene = ();
my %gene_names_existing_in_ensembl = (); my %ens_transcript_id_to_gene_name_lookup = ();
my $transcript_id; my $biotype; my $gene_id; my $gene_name; my $desc;
open(IN,$in_file1) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $transcript_id = $1; $gene_id = ''; $gene_name = '';
		  if ($transcript_id =~ /^(.*?)\.\d+$/)
			{ $transcript_id = $1; }
		  if ($line =~ /^.*? gene\_biotype\:(.*?) transcript\_biotype.*?$/)
			{ $biotype = $1; }
		  if ($line =~ /^.*? gene\:(.*?)\.\d+ gene\_biotype\:.*?$/)
			{ $gene_id = $1; }
		  if (($line =~ /^.*? gene\_symbol\:(.*?) description.*?$/) || ($line =~ /^.*? gene\_symbol\:(.*?)$/))
			{ $gene_name = $1; }
		  if ($line =~ /^.*? description\:(.*?) \[Source.*?$/)
			{ $desc = $1; }
		  if ($gene_name eq '') { $gene_name = $gene_id; }
		  $gene_names_existing_in_ensembl{$gene_name}++;
		}
	  else
		{ if ($biotype eq 'protein_coding')
			{ $ensembl_seqs{$transcript_id} .= uc($line);
			  $ensembl_seqs_per_gene{$gene_name}{$transcript_id} .= uc($line);
			}
		  $ens_transcript_id_to_gene_name_lookup{$transcript_id}{gene_name} = $gene_name;
		  $ens_transcript_id_to_gene_name_lookup{$transcript_id}{gene_id}   = $gene_id;
		  $ens_transcript_id_to_gene_name_lookup{$transcript_id}{desc}      = $desc;
		}
	}
close(IN) or die $!;

# STORE PROTEIN-CODING TRANSCRIPTS (NCBI REFSEQS)
my %ncbi_seqs = (); my %usable_ncbi_transcripts = (); my %ncbi_transcript_id_details = ();
open(IN,$in_file2) or die $!;
while(<IN>)
	{ my $line = $_; chomp($line);
	  if ($line =~ /^\>(.*?) .*?$/)
		{ $transcript_id = $1;
		  if ($transcript_id =~ /^(.*?)\.\d+$/)
			{ $transcript_id = $1;
			  if ($line =~ /^\>.*? $species (.*) \((.*?)\)\, (.*?)mRNA$/) # FILTER: exclude non-mRNA genes (such as 'misc_RNA' or 'partial mRNA')
				{ my $desc = $1; my $gene_name = $2; my $type = $3;
				  if ($gene_name =~ /^C(\d+)ORF(\d+)$/) { $gene_name = "C"."$1"."orf"."$2"; }
				  if (($type !~ /partial/) && ($type !~ /misc\_/)) # this qualifier is necessary so that we can retain sequences with this type of header: ">XM_015294055.1 PREDICTED: Gallus gallus CD163 molecule (CD163), transcript variant X1, mRNA"
					{ next if (!(exists($ncbi_transcript_id_to_gene_id_lookup{$transcript_id})));
					  my $gene_id = $ncbi_transcript_id_to_gene_id_lookup{$transcript_id};
					  $ncbi_transcript_id_details{$transcript_id}{gene_name} = $gene_name;
					  $ncbi_transcript_id_details{$transcript_id}{gene_id} = $gene_id;
					  $ncbi_transcript_id_details{$transcript_id}{desc} = $desc;
					  $usable_ncbi_transcripts{$transcript_id}++;
					}
				}
			}
		}
	  else
		{ $ncbi_seqs{$transcript_id} .= uc($line) unless (!(exists($usable_ncbi_transcripts{$transcript_id})));
		}
	}
close(IN) or die $!;

# WHICH TRANSCRIPT SEQUENCES ARE IDENTICAL IN BOTH THE ENSEMBL AND REFSEQ DATASETS? WE'LL IDENTIFY THESE SO THAT WE MAY LATER OUTPUT IN OUR FINAL REFERENCE TRANSCRIPTOME ONLY THE ENSEMBL GENES, SUPPLEMENTING IT WITH THOSE REFSEQS SUPPORTING GENE NAMES NOT YET PRESENT IN ENSEMBL
my %all_seqs = ();
while((my $transcript_id,my $seq)=each(%ensembl_seqs)) { $all_seqs{$seq}{ens}  = $transcript_id; }
while((my $transcript_id,my $seq)=each(%ncbi_seqs))    { $all_seqs{$seq}{ncbi} = $transcript_id; }
my %ncbi_refseqs_to_retain = ();
while((my $seq,my $irrel)=each(%all_seqs))
	{ if ((exists($all_seqs{$seq}{ncbi})) and (!(exists($all_seqs{$seq}{ens}))))
		{ my $ncbi_refseq_id = $all_seqs{$seq}{ncbi};
		  my $gene_name = $ncbi_transcript_id_details{$ncbi_refseq_id}{gene_name};
		  if (exists($gene_names_existing_in_ensembl{$gene_name}))
			{ # FILTER: at this point, an NCBI mRNA sequence for this gene exists, whereas a corresponding Ensembl CDS does not - but the mRNA can include UTRs. As such, if any CDS is encapsulated by this mRNA, it will only encode the same protein as that CDS... and so it should not be retained.
			  my $this_mrna_includes_an_ensembl_cds = 0;
			  while((my $ens_transcript_id,my $cds)=each(%{$ensembl_seqs_per_gene{$gene_name}}))
				{ if ($seq =~ /^.*?$cds.*?$/)
					{ $this_mrna_includes_an_ensembl_cds++;
					}
				}
			  if ($this_mrna_includes_an_ensembl_cds == 0)
				{ $ncbi_refseqs_to_retain{$ncbi_refseq_id}++ unless ($ncbi_refseq_id =~ /^NM\_(\d+)$/);
				}
			}
		  else
			{ $ncbi_refseqs_to_retain{$ncbi_refseq_id}++ unless ($ncbi_refseq_id =~ /^NM\_(\d+)$/); } # FILTER: we will include in the final reference transcriptome RefSeq IDs only if they (a) can be assigned to a gene not yet named in Ensembl (e.g. CD163), or (b) are a novel mRNA for an existing gene in Ensembl
		}
	}

# CREATE REFERENCE TRANSCRIPTOME
my %ens_genes_seen = (); my %ncbi_genes_seen = (); my %tx2gene = ();
my $ens_transcripts_seen = 0; my $ncbi_transcripts_seen = 0;
my @transcript_ids = ();
while((my $transcript_id,my $seq)=each(%ensembl_seqs))
	{ push(@transcript_ids,$transcript_id); }
my @sorted_transcript_ids = sort {$a cmp $b } @transcript_ids;
for(my $x=0;$x<@sorted_transcript_ids;$x++)
	{ my $pc = sprintf("%.2f",(($x/$#sorted_transcript_ids)*100));
	  print "parsing Ensembl transcripts: $x of $#sorted_transcript_ids ($pc%)\n";
	  my $transcript_id = $sorted_transcript_ids[$x];
	  my $seq = $ensembl_seqs{$transcript_id};
	  next if (!(exists($ens_transcript_id_to_gene_name_lookup{$transcript_id})));
	  my $length = length($seq);
	  print OUT ">$transcript_id\n$seq\n";
	  $ens_transcripts_seen++;
	  my $desc      = $ens_transcript_id_to_gene_name_lookup{$transcript_id}{desc};
	  my $gene_id   = $ens_transcript_id_to_gene_name_lookup{$transcript_id}{gene_id};
	  my $gene_name = $ens_transcript_id_to_gene_name_lookup{$transcript_id}{gene_name};
	  $ens_genes_seen{$gene_name}++;
	  my $gene_loc = '';
	  if    (exists($ens_gene_locs{$gene_id}))  { $gene_loc =  $ens_gene_locs{$gene_id}; }
	  elsif (exists($ncbi_gene_locs{$gene_id})) { $gene_loc = $ncbi_gene_locs{$gene_id}; }
	  else
		{ print "ERROR: unable to find a location for gene $gene_id\n"; exit 1; }
	  push(@{$tx2gene{$gene_name}{$transcript_id}},"$gene_name\t$gene_id\t$transcript_id\t$desc\t$gene_loc");
	}
my %non_ens_genes_seen_in_ncbi = (); my %ens_genes_with_novel_transcripts = ();
@transcript_ids = ();
while((my $transcript_id,my $seq)=each(%ncbi_seqs))
	{ push(@transcript_ids,$transcript_id); }
@sorted_transcript_ids = sort {$a cmp $b } @transcript_ids;
for(my $x=0;$x<@sorted_transcript_ids;$x++)
	{ my $pc = sprintf("%.2f",(($x/$#sorted_transcript_ids)*100));
	  print "parsing NCBI transcripts: $x of $#sorted_transcript_ids ($pc%)\n";
	  my $ncbi_refseq_id = $sorted_transcript_ids[$x];
	  if (exists($ncbi_refseqs_to_retain{$ncbi_refseq_id}))
		{ my $seq = $ncbi_seqs{$ncbi_refseq_id};
	      # remove UTRs from mRNA by identifying the longest ORF and excising everything outside it
		  my %seq_by_frame = ();
		  for(my $frame=1;$frame<=6;$frame++)
			{ my $dna_seq = $seq;
			  my $translate_frame_2 = 0; my $translate_frame_3 = 0;
			  my $comp_seq = 0;
			  if (($frame == 4) || ($frame == 5) || ($frame == 6)) { $comp_seq = 1; }
			  if (($frame == 2) || ($frame == 5))
				{ $dna_seq =~ s/^\w{1}// unless ($translate_frame_2 > 0);
				  $translate_frame_2++;
				}
			  if (($frame == 3) || ($frame == 6))
				{ $dna_seq =~ s/^\w{2}// unless ($translate_frame_3 > 0);
				  $translate_frame_3++;
				}
			  if ($comp_seq == 1) { $dna_seq =~ tr/[ATCG]/[TAGC]/; }
			  $seq_by_frame{$frame} = $dna_seq;
			}
		  my @possible_orfs = ();
		  my %frames = qw {1 +1 2 +2 3 +3 4 -1 5 -2 6 -3};
		  my $last_end = 0; my $last_orf_length_seen = 0;
		  for(my $frame=1;$frame<=6;$frame++)
			{ my $exact_frame = $frames{$frame};
			  my $this_seq = $seq_by_frame{$frame};
			  my @seq = split(//,$this_seq);
			  my @orf_starts = (); my @orf_ends = ();
			  my $codon_num = 1;
			  for(my $x=0;$x<@seq;$x++)
				{ my $nt1 = $seq[$x];
				  next if ( (!(defined($seq[$x+1]))) || (!(defined($seq[$x+2]))) );
				  my $nt2 = $seq[$x+1]; my $nt3 = $seq[$x+2];
				  my $codon = "$nt1"."$nt2"."$nt3";
				  if ($codon eq 'ATG')
					{ my $start = $x+1;
					  if 	(($frame == 2) || ($frame == 5)) { $start = $start+1; }
					  elsif (($frame == 3) || ($frame == 6)) { $start = $start+2; }
					  push(@orf_starts,$start);
					}
				  elsif (($codon eq 'TAA') or ($codon eq 'TAG') or ($codon eq 'TGA'))
					{ my $end = $x+3;
					  if 	(($frame == 2) || ($frame == 5)) { $end = $end+1; }
					  elsif (($frame == 3) || ($frame == 6)) { $end = $end+2; }
					  push(@orf_ends,$end);
					}
				  $codon_num++;
				  $last_end = $x+3;
				  $x = $x+2;
				}
			  push(@orf_ends,$last_end);
			  next if ($#orf_starts == -1);
			  
			  # what is the longest distance between any two values of the ORF start and ORF end arrays? Viable ORFs are those with $orf_end nearest to $orf_start. We will later consider only the longest viable ORF.
			  my @orfs = ();
			  foreach my $orf_start (@orf_starts)
				{ my @nearest_orf_end = ();
				  foreach my $orf_end (@orf_ends)
					{ my $orf_length = ($orf_end-$orf_start)+1;
					  push(@nearest_orf_end,[$orf_length,$orf_end]) unless ($orf_length < 1);
					}
				  my @sorted_nearest_orf_end = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, $_->[0]] } @nearest_orf_end;
				  my $orf_length = $sorted_nearest_orf_end[0][0]; my $orf_end = $sorted_nearest_orf_end[0][1];
				  push(@orfs,[$orf_length,$orf_start,$orf_end]); # this is a complete ORF (longest distance between any start and the end)
				}
			  my @sorted_orfs = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @orfs;
			  my $orf_start = $sorted_orfs[0][1]; my $orf_end = $sorted_orfs[0][2];
			  my $orf_length = ($orf_end-$orf_start)+1;
			  my $orf = substr($seq,$orf_start-1,$orf_length);
			  my $strand = '';
			  if ($exact_frame < 0) { $strand = '-'; } else { $strand = '+'; }
			  if ($exact_frame < 0) { $orf =~ tr/[ATCGN]/[TAGCN]/; }
			  $exact_frame =~ s/[\-\+]//g;
			  my @orf = split(//,$orf);
			  my $translated_orf = '';
			  for(my $x=0;$x<@orf;$x++)
				{ my $nt1 = $orf[$x];
				  next if ( (!(defined($orf[$x+1]))) || (!(defined($orf[$x+2]))) );
				  my $nt2 = $orf[$x+1]; my $nt3 = $orf[$x+2];
				  if ($nt1 !~ /^[ATCGN]$/) { $nt1 = 'N'; } # replace all ambiguity characters (R, Y, S, W, K, M, B, D, H, V) with N
				  if ($nt2 !~ /^[ATCGN]$/) { $nt2 = 'N'; }
				  if ($nt3 !~ /^[ATCGN]$/) { $nt3 = 'N'; }
				  my $codon = "$nt1"."$nt2"."$nt3";
				  my $aa = $trans{$codon};
				  $translated_orf .= $aa;
				  $x = $x+2;
				}
			  push(@possible_orfs,[$orf_length,$orf_start,$orf_end,$orf,$translated_orf]) if ($orf_length >= $last_orf_length_seen);
			  $last_orf_length_seen = $orf_length if ($orf_length > $last_orf_length_seen);
			}
		  my @sorted_possible_orfs = map { $_->[0] } sort { $b->[1] <=> $a->[1] } map { [$_, $_->[0]] } @possible_orfs;
		  my $orf_seq = $sorted_possible_orfs[0][3];
		  my $aa_seq  = $sorted_possible_orfs[0][4];
		  my $length  = length($orf_seq);
		  print OUT ">$ncbi_refseq_id\n$orf_seq\n";
		  $ncbi_transcripts_seen++;
		  my $gene_name = $ncbi_transcript_id_details{$ncbi_refseq_id}{gene_name};
		  my $gene_id   = $ncbi_transcript_id_details{$ncbi_refseq_id}{gene_id};
		  my $desc      = $ncbi_transcript_id_details{$ncbi_refseq_id}{desc};
		  $ncbi_genes_seen{$gene_name}++;
		  if (exists($ens_genes_seen{$gene_name}))
			{ $ens_genes_with_novel_transcripts{$gene_name}++; }
		  else
			{ $non_ens_genes_seen_in_ncbi{$gene_name}++; }
		  my $gene_loc = '';
		  if    (exists($ens_gene_locs{$gene_id}))  { $gene_loc =  $ens_gene_locs{$gene_id}; }
		  elsif (exists($ncbi_gene_locs{$gene_id})) { $gene_loc = $ncbi_gene_locs{$gene_id}; }
		  else
			{ print "ERROR: unable to find a location for gene $gene_id\n"; exit 1; }
		  push(@{$tx2gene{$gene_name}{$ncbi_refseq_id}},"$gene_name\t$gene_id\t$ncbi_refseq_id\t$desc\t$gene_loc");
		}
	}
close(OUT) or die $!;

my @gene_names = ();
while((my $gene_name,my $irrel)=each(%tx2gene))
	{ push(@gene_names,$gene_name); }
my @sorted_gene_names = sort {$a cmp $b} @gene_names;
foreach my $gene_name (@sorted_gene_names)
	{ my @transcript_ids = ();
	  while((my $transcript_id,my $irrel)=each(%{$tx2gene{$gene_name}}))
		{ push(@transcript_ids,$transcript_id); }
	  my @sorted_transcript_ids = sort {$a cmp $b} @transcript_ids;
	  foreach my $transcript_id (@sorted_transcript_ids)
		{ foreach my $line (@{$tx2gene{$gene_name}{$transcript_id}})
			{ print TX2GENE "$line\n"; }
		}
	}
close(TX2GENE) or die $!;
my $ens_genes_seen  = scalar keys %ens_genes_seen;
my $ncbi_genes_seen = scalar keys %ncbi_genes_seen;
my $non_ens_genes_seen_in_ncbi = scalar keys %non_ens_genes_seen_in_ncbi;
my $ens_genes_with_novel_transcripts = scalar keys %ens_genes_with_novel_transcripts;
print "reference transcriptome contains $ens_transcripts_seen Ensembl transcripts ($ens_genes_seen genes) and $ncbi_transcripts_seen NCBI transcripts ($ncbi_genes_seen genes)\n";
print "of the $ncbi_genes_seen NCBI genes, $ens_genes_with_novel_transcripts are Ensembl genes with novel transcripts and $non_ens_genes_seen_in_ncbi are specific to NCBI\n";

exit 1;