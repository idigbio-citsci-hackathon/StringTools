#!/usr/bin/perl -w
use strict;

if (!@ARGV) {
  warn "usage:\n $0 'some names'\n";
  while (<DATA>) {
    chomp;
    my @list = split /\s*\|\s*/;
    my $in = shift @list;
    my @out = list_split($in);
    
    my $expected = join ' | ',@list;
    my $actual = join '| ',@out;
    if ($expected ne $actual) {
      warn "FAIL: expected '$expected'; actual '$actual'\n";
    }
  }
} else {
  for my $name (@ARGV) {
    my @split = list_split($name);
    print join "\t",@split;
    print "\n";
  }
}
sub list_split {
  my $in = shift;
  return ($in);
}

__DATA__
Stuttgart | Stuttgart
James D. Ray Jr., C. Earle Smith, Jr., Olga Lakela, Jackie Patman | James D. Ray Jr. | C. Earle Smith, Jr. | Olga Lakela | Jackie Patman