#!/usr/bin/env perl
#


$argc = @ARGV;
if ($argc == 0) {
    print "usage: cat aln1.sam | sam_add_rg.pl <read group id> <sample name>\n";
    print "changes all alignments to have RG tags matching the new sample name.\n";
    print "prints SAM format out stdout\n";
    exit(0);
}

$ID = $ARGV[0];
$SM = $ARGV[1];

$in_header = 1;

while (<STDIN>) {
    if ($_ =~ /^@/) {
        print $_;
        next;
    } else {
        # add the new RG group to the end of the header
        if ($in_header) {
            print "\@RG\tID:$ID\tSM:$SM\n";
            $in_header = 0;
        }
        if ($_ =~ /\tRG:Z:(.+?)/) {
            $_ =~ s/\tRG:Z:.+?([\t\n])/\tRG:Z:$ID$1/;
        } else {
            $_ =~ s/\n/\tRG:Z:$ID\n/;
        }
        print $_;
    }
}

