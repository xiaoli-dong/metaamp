#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper ;
use JSON;

my $filename = $ARGV[0];
my $json_text;
{
  local $/;
  my $fh;
  open $fh, $filename;
  $json_text = <$fh>;
  close $fh;
}
#my @json = @{decode_json($json_text)};
my @json = decode_json($json_text);
print Dumper(@json);
