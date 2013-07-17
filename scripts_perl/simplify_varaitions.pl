#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use File::Basename;

my $help;
my $candidates;

GetOptions(
	'candidates=s'	=> \$candidates,
	'help|?' 		=> \$help
) or pod2usage(2);


system("grep -v random $candidates | grep -v chrUn > tmp_noRandomChr_noChrUn.bed");

open(FILE, "tmp_noRandomChr_noChrUn.bed");



my $line = <FILE>;
chomp $line;
my @variation = (split /\t/, $line);

my $current_ID_fc    = $variation[0];
my $current_start_fc = $variation[1];
my $current_end_fc   = $variation[2];

my $current_ID_sc    = $variation[3];
my $current_start_sc = $variation[4];
my $current_end_sc   = $variation[5];

my @CHRspec_Variations = ();
push @CHRspec_Variations, $line;

while($line = <FILE>) {
    chomp $line;
    @variation = (split /\t/, $line);
    my $new_ID_fc    = $variation[0];
    my $new_start_fc = $variation[1];
    my $new_end_fc   = $variation[2];
    
    my $new_ID_sc    = $variation[3];
    my $new_start_sc = $variation[4];
    my $new_end_sc   = $variation[5];
    
    if($new_ID_fc ne $current_ID_fc ) {
        #### process the chr and look for intervals that are close
        my @NewVariations = ();
        if(@CHRspec_Variations > 1) {
            ## check if we can merge some lines (when multiple links are present on the same pair)
            my %secondChrHits;
            foreach my $entry (@CHRspec_Variations) {
                my @entry = (split /\t/, $entry);
                push @{$secondChrHits{$entry[3]}}, $entry;
            }
            
            print "Hits on $new_ID_fc\n";
            foreach my $secondHit (keys %secondChrHits) {
                my @vectorHits = @{$secondChrHits{$secondHit}};
                
                foreach my $hit (@vectorHits) {
	                print "\t$hit\n";
                }
                print "#######\n";
                
            }
            
        } else {
            print $CHRspec_Variations[0]."\n";
        }
        
        
        @CHRspec_Variations = ();
        push @CHRspec_Variations, $line;
       
        $current_ID_fc    = $new_ID_fc;
        $current_start_fc = $new_start_fc;
        $current_end_fc   = $new_end_fc;
        
        $current_ID_sc    = $new_ID_sc;
        $current_start_sc = $new_start_sc;
        $current_end_sc   = $new_end_sc;
        
        
        
    } else {
         push @CHRspec_Variations, $line;
    }
}

foreach my $entry (@CHRspec_Variations) {
    print $entry."\n";
}



