use strict;
use warnings;
use JSON qw(decode_json);

# REQUIREMENTS
my $json = $ARGV[0];

open(JSON,'<',$json) or die $!;
local $/;
my $fastp = decode_json(<JSON>);
close(JSON);
my $avg_read_length = $fastp->{'summary'}{'before_filtering'}{'read1_mean_length'};
print "$avg_read_length";