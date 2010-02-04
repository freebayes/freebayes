#!/usr/bin/perl

use strict;

my $n = shift;
my $th = shift;

my $sum = 0;
for(my $i = 1; $i < $n; $i++) {
    $sum += 1/$i;
}

print "SNP rate=", $th * $sum, " or 1 SNP per ", 1 / ($sum*$th), " bases\n";
