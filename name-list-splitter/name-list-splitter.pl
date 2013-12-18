#!/usr/bin/perl -w
use strict;

### CLI and test runner:

if (!@ARGV) {
  warn "usage:\n $0 'some names'\n";
  my $fail = 0;
  while (<DATA>) {
    chomp;
    s/#.*$//;
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

### Raisin d'etre:

sub list_split {
  my $in = shift;
  $in =~s/^collected by[^\w]*//i;
  my @tokens = grep {$_} # drop empty tokens
    split m{
      (?<=[,;]) | # keep comma with preceding.
      (?:\s+(?:with|and|&)\s+) |
      (?:\s+w/)
    }x, $in; 
  my @out;
  while (@tokens) {
    my $token = shift @tokens;
    if ($token =~ /^\s*jr\.?\s*,?\s*$/i && @out) {
      $out[scalar @out - 1] .= $token;
    } else {
      push @out, $token;
    }
  }
  @out = map {s/[,;]$//;$_} @out;
  return @out;
}

### Test data:

# Note that this is just splits and maps names: no other clean-up is in scope.
# "det" or "phd" are not touched.

__DATA__
Stuttgart | Stuttgart
James D. Ray Jr. | James D. Ray Jr.
C. Earle Smith, Jr. | C. Earle Smith, Jr.
Olga Lakela, Jackie Patman | Olga Lakela | Jackie Patman
James D. Ray Jr., C. Earle Smith, Jr., Olga Lakela, Jackie Patman | James D. Ray Jr. | C. Earle Smith, Jr. | Olga Lakela | Jackie Patman
D. S. Correll and Helen B. Correll | D. S. Correll | Helen B. Correll
R. K Godfrey with Angus Gholson | R. K Godfrey | Angus Gholson
H. Maurushat; V. Sullivan, C. Hudson; D. Wise, R.K. Godfrey | H. Maurushat | V. Sullivan | C. Hudson | D. Wise | R.K. Godfrey
Collected by: R. K. Godfrey | R. K. Godfrey
william p adams | william p adams
A. Gholson Jr w/Wilson Baker | A. Gholson Jr | Wilson Baker
D. B. Ward, with H. F. Decker | D. B. Ward | H. F. Decker
R. K. Godfrey (det.) & Richard D. Houk | R. K. Godfrey (det.) | Richard D. Houk

### These don't pass yet:

# Cecil R. Slaughter, Ph.D. | Cecil R. Slaughter, Ph.D.
# R, Kral & P.L. Redfearn | R, Kral | P.L. Redfearn
# D. B. & S. S. Ward | D. B. Ward | S. S. Ward
# R> K> Godfrey with Robt. & John Lazor | R> K> Godfrey | Robt. Lazor | John Lazor
# R, K, Godfrey | R, K, Godfrey # This was a real transcription!
# R. K. Godfrrey with Robt. & John Lazor | R. K. Godfrrey | Robt. Lazor | John Lazor
