#!/usr/bin/perl -w
use strict;


sub list_split {
  my $in = shift;
  my @tokens = split /(?<=,)/, $in; # keep comma with preceding.
  my @out;
  while (@tokens) {
    my $token = shift @tokens;
    if ($token =~ /^\s*jr\.?\s*,?\s*$/i && @out) {
      $out[scalar @out - 1] .= $token;
    } else {
      push @out, $token;
    }
  }
  @out = map {s/,$//;$_} @out;
  return @out;
}


if (!@ARGV) {
  warn "usage:\n $0 'some names'\n";
  my $fail = 0;
  while (<DATA>) {
    chomp;
    next if /^#/;
    next unless /./;
    my @list = split /\s*\|\s*/;
    my $in = shift @list;
    my @out = list_split($in);
    
    my $expected = join ' | ',@list;
    my $actual = join ' | ',@out;
    $expected =~s/\s+\|\s+/ | /g;
    $actual =~s/\s+\|\s+/ | /g;
    if ($expected ne $actual) {
      warn "FAIL: expected '$expected'; actual '$actual'\n";
      $fail++;
    }
  }
  if (!$fail) {
    warn "(all tests pass.)\n";
  }
} else {
  for my $name (@ARGV) {
    my @split = list_split($name);
    print join "\t",@split;
    print "\n";
  }
}


__DATA__
Stuttgart | Stuttgart
James D. Ray Jr. | James D. Ray Jr.
C. Earle Smith, Jr. | C. Earle Smith, Jr.
Olga Lakela, Jackie Patman | Olga Lakela | Jackie Patman
James D. Ray Jr., C. Earle Smith, Jr., Olga Lakela, Jackie Patman | James D. Ray Jr. | C. Earle Smith, Jr. | Olga Lakela | Jackie Patman