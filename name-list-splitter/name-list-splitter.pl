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
  $in =~s/(?<![A-Z])([A-Z]),/$1./g; # There were lots of commas after initials.
  $in =~s/\.+/./g;
  my @tokens = grep {$_} # drop empty tokens
    split m{
      (?<=[,;]) | # keep comma with preceding.
      (?:\s+\W?(?:with|and)\W?\s+) |
      (?:\s+w/) |
      (?:\+|&)
    }xi, $in; # I don't have positive examples where case-insensitive is actually necessary.
  my @out;
  while (@tokens) {
    my $token = shift @tokens;
    if ($token =~ /
        ^\s*(
          jr\.? |
          sr\.? |
          ph\.?d\.?
        )\s*,?\s*$/xi && @out) {
      $out[scalar @out - 1] .= $token;
    } else {
      push @out, $token;
    }
  }
  @out = map {s/[,;]$//;$_} @out;
  return @out;
}

### Test data:

# Note that this is just splits name strings, and associates stray initials with the right last name.
# No other clean-up is in scope.
# "det" or "phd" are not touched.

# No pretensions of making this a general solution: in particular
#   -- we do want to handle commas-that-should-have-been-periods, since we do see that in these transcriptions
#   -- and we don't worry about last-name-first style, since I haven't seen any of those. 

# All test cases are taken from real crowd-sourced data.

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
Cecil R. Slaughter, Ph.D. | Cecil R. Slaughter, Ph.D.
George Eiten, Liene T. Eiten, Gil M. Felippe & J.M. de Freitas Campes | George Eiten | Liene T. Eiten | Gil M. Felippe | J.M. de Freitas Campes
Tim Reeves, Sr | Tim Reeves, Sr
John W/ Thieret | John | Thieret
A. Gholson, Jr., w/ Dr. Bob Godfrey & D. C. Vickers | A. Gholson, Jr. | Dr. Bob Godfrey | D. C. Vickers
Gholson, Jr., with Dr Bob Godfrey and D.C. Vickens | Gholson, Jr. | Dr Bob Godfrey | D.C. Vickens
Charles T. Bryson (and) Will McDearman | Charles T. Bryson | Will McDearman
S Leonard + H McAninch | S Leonard | H McAninch
Mrs. C. W. Hanes | Mrs. C. W. Hanes
Coll. Ame Garthwright | Coll. Ame Garthwright
W.H. Wagner , Jr. | W.H. Wagner , Jr.

### Comma which should be period:

R, Kral & P.L. Redfearn | R. Kral | P.L. Redfearn
R, K, Godfrey | R. K. Godfrey
W. A, Sliveus | W. A. Sliveus
R,K, Godfrey & J.P, Gillespie | R.K. Godfrey | J.P. Gillespie
R. K,. Godfrey and Richard D. Houk | R. K. Godfrey | Richard D. Houk

### ACK!! No punctuation between names:

# M.B. H.L. | M.B. H.L.
# D.R.Windler B.r. Sinor | D.R.Windler | B.r. Sinor
# S.W.Leonard D. Culwell M.Ripperton | S.W.Leonard | D. Culwell | M.Ripperton
# R. K. Godfrey Richard D. Houk | R. K. Godfrey | Richard D. Houk

### Name distribution:

# Nancy Craft Coile, w/ Robert, Danielle & Robbie Coile | Nancy Craft Coile | Robert Coile | Danielle Coile | Robbie Coile
# P.J. Crutchfeld & Laura & Thomas Crutchfield | P.J. Crutchfeld | Laura Crutchfield | Thomas Crutchfield
# R. K. Godfrrey with Robt. & John Lazor | R. K. Godfrrey | Robt. Lazor | John Lazor
# D. B. & S. S. Ward | D. B. Ward | S. S. Ward

### Slashes:

# Michel G. Lelong / Ken Rogers | Michel G. Lelong | Ken Rogers
# SW Leonard/D. Culwell/M. Ripperton | SW Leonard | D. Culwell | M. Ripperton
# A. H. Curtiss/ det. C. B. Heiser, Jr.
# R/K/ Godfrey & John Morrill
# Robert & Mabel Kral/

### Parens:

# John Mayberg (By W)
# (Mary L. Leigh) J. Rowntrey
# A.R. Diamond (w. J.D. Freeman)
# David Hall (w/ Gary Schultz)
# R K Godfrey ( Shirley Mah Kooyman 1980)
# (Karl, Godfrey 1958); R K Godfrey 1976
# Loran C Anderson ( Scott Sundberd 1987)
# (MARY L. LEIGH) J. ROUNTREY

### Dashes:

# H.H. Iltis - D. Parker
# ROBERT K GODFREY - R W SIMONS - ANGUS GHOLSON
# Robert F- Thorne
# Steve L. Orzell- Edwin L. Bridges
# A GHOLSON, JR - SUSANNE COOPER - WILSON BAKER
# H Maurushat - V. Sullivan, C. Hudson
# REGINALD TOSE-INNES & BARTON H. WARNOCK
# R.K. Godfrey w/- Christopher Campbell
# Barton H. warnock, Reginald Rose-Innes
# Marie-Victorin, Rolland-Germain, Marcel Raymond

### Other:

# R> K> Godfrey with Robt. & John Lazor | R> K> Godfrey | Robt. Lazor | John Lazor
# Lytton J. Musselman / Elizabeth R. Musselman
# Judith Canne and Jose Schunkell (or Schunkey?)
# ?Rufus Crane
# A.H.S.F.
# Lloyd T. (?Y.?) Card (?Cart?)
# R:D. Houk ans R:K: Godfrey
# Bruce Hansen with T.&B. Cochrane, C.S. Keller & M. Waterway
# Loran C: Anderson

### ???:

# Loran C. Anderson w/Gil Nelson R>K> Godfrey Herbarium (FSU)
# Ann F-Johnson
# D.B. Ward 3-19, with BTY 421
# Betts, Thealcald Jonas, Baker.
# Donald Eves (1956) and L. J. Ultal (1981)
# Lafon & Gray Bill

