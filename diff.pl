#!/usr/bin/env perl

use strict;
use warnings;

sub readFile
{
  my ($file, $h) = @_;
  open INPUT, $file or die qq|Failed to open $file, died|;
  while (<INPUT>) {
    chop;
    my ($run, $lumi, $event, $mass4l, $massz1, $massz2) = (split /:/)[0,1,2,3,4,5]; 
    $h->{$run.":".$lumi.":".$event.":".$mass4l.":".$massz1.":".$massz2} = $_;
  }
  close INPUT;
}
sub _main {
  my ($infile1, $infile2) = @_;
  my $dict1 = {};
  readFile($infile1, $dict1);

  my $dict2 = {};
  readFile($infile2, $dict2);

  print join("", ">>> Lines present in <", $infile1, "> not in <", $infile2, ">"), "\n";
  for my $key (sort keys %$dict1) {
    print $dict1->{$key}, "\n" unless exists $dict2->{$key};
  }

  print join("", ">>> Lines present in <", $infile2, "> not in <", $infile1, ">"), "\n";
  for my $key (sort keys %$dict2) {
    print $dict2->{$key}, "\n" unless exists $dict1->{$key};
  }
}
_main($ARGV[0], $ARGV[1]);
__END__
