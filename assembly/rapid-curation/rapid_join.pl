#!/software/bin/perl -w
# kj2 08.10.2020
# jmdw added in hap flag and changed output suffix for files

use strict;
use Bio::SeqIO;
use Getopt::Long;
no warnings 'uninitialized'; # turn that off for debugging

my $tpf;
my $fa;
my $out;
my $csvin;
my $hap;
my $help;

GetOptions (
    "fa:s"  => \$fa,
    "tpf:s" => \$tpf,
    "out:s" => \$out,
    "csv:s" => \$csvin,
    "hap"   => \$hap,
    "h"     => \$help,
    "help"  => \$help,
);

if (($help) || (!$fa) || (!$tpf)) {
    print "This script takes an original assembly fasta file, a one-line per chromosome pre-csv file and a TPF file generated with split.pl and creates the assembly from the TPF file\n";
    print "Usage:\n";
    print "perl join.pl -fa <fasta>\n";
    print "             -tpf <tpf>\n";
    print "             -csv <pre-csv>\n";
    print "             -out <outfile_fasta_prefix> \# unless specified the output will be written to <fasta>.curated.fasta\n";
    print "             -hap \# optional use only if generating haplotigs fasta\n";
    print "             -h/help # this message\n";   
    exit(0);
}

my $newout;
if ($out) {
    $newout = $out;
}
else {
    ($newout) = ($fa =~ /(\S+)\.fa/);
    $out = ${newout}.".intermediate.fa" unless ($out);
}

my $seqout;
if ($hap) {
	$seqout = Bio::SeqIO->new(-format => 'fasta',
				-file => "> ${newout}.additional_haplotigs.unscrubbed.fa");
}
else {
	$seqout = Bio::SeqIO->new(-format => 'fasta',
				-file => "> ${newout}.inter.fa");
}

# storing sequences
my %seqhash;
my $seqin  = Bio::SeqIO->new('-format' => 'fasta',
                             '-file'   => $fa);
while (my $seqobj = $seqin->next_seq()) {
    die("ERROR: you are trying to run rapid_join.pl on a full curation file\n") if $seqobj->display_id =~/_ctg1/;
    $seqhash{$seqobj->display_id} = $seqobj;
}

open(TPF,$tpf);
my $seqobj;
my $lastscaff;
my $seqstring;
my %seen;
my $no;
my $lastgap;
my %created;
my %assembly;

while (<TPF>) {
    $no++;
    my $line = $_;
    next if $line =~/^\n/;
    chomp;
    
    # gaps
    if (/^GAP/) {
        # no gaps at start of component
        die("Sequence can't begin with gap in line $no\n") unless ($lastscaff);

        my $length;
        if (/^GAP\s+\S+\s+(\d+)/) {
            $length = $1;
        }
        else {
            $length = 200;
        }
        $seqstring .= "N"x$length;
        $lastgap++;
    }
    # sequence
    else {
        my ($undef,$ctg,$scaff,$ori) = split;
        my ($ctgname,$start,$end) = ($ctg =~ /(\S+)\:(\d+)-(\d+)/);
        
        # check whether scaffold name was used previously
        die("ERROR: attempt to reuse scaffold name $scaff, line $no\n") if (exists $created{$scaff});
        
        # check whether sequence is being reused
        push @{$seen{$ctgname}}, [$start,$end];
        
        # continuation of scaffold
        if ($lastscaff &&($scaff eq $lastscaff)) {
            
            undef $lastgap;
            my $oldobj = $seqhash{$ctgname};
            my $subseq;
            eval {
                $subseq = $oldobj->subseq($start,$end);
            };   
            die("ERROR: $scaff:$start-$end does not exist\n") if ($@);
            
                 
            # revcom
            if ($ori eq "MINUS") {
                $subseq = reverse($subseq);
                $subseq =~ tr/atcgATCG/tagcTAGC/;
            }
            
            $seqstring .= $subseq;
        }
        
        # new scaff
        else {
            
            # check whether there's a gap at start of last scaffold
            die("ERROR: Previous component ended in gap, line ",$no-1,"\n") if ($lastgap);
            
            # out with last
            if ($lastscaff) {
                $seqobj->seq($seqstring);
                $seqobj->display_id($lastscaff);
                $seqout->write_seq($seqobj);
                $assembly{$seqobj->display_id($lastscaff)} = $seqobj;
                $created{$lastscaff}= length($seqstring);
                #print "Created $lastscaff with ", length($seqstring),"\n";

                undef $seqobj;
                undef $seqstring;
            }
            
            $seqobj = Bio::Seq->new();
            my $oldobj = $seqhash{$ctgname};
            $seqobj->display_id($scaff);
            
            my $subseq = $oldobj->subseq($start,$end);
            
            # revcom
            if ($ori eq "MINUS") {
                $subseq = reverse($subseq);
                $subseq =~ tr/atcgATCG/tagcTAGC/;
            }
            
            $seqstring .= $subseq;
            $lastscaff = $scaff;
            
            undef $lastgap;
        
        }
    }
}
# out with last
$seqobj->seq($seqstring);
$seqobj->display_id($lastscaff);
$seqout->write_seq($seqobj);
$assembly{$seqobj->display_id($lastscaff)} = $seqobj;
$created{$lastscaff}= length($seqstring);
#print "Created $lastscaff with ", length($seqstring),"\n";

# check for overlaps
foreach my $ctg (keys %seen) {
    my @sorted = sort {$a->[0] <=> $b->[0]} @{$seen{$ctg}};
    my $max;
    foreach my $aref (@sorted) {
        my $start = $aref->[0];
        my $end = $aref->[1];
        my $tick = 0;
        if ($max) {
            if ($start < $max) {
                $tick++;    
            }
            elsif (($start < $max) && ($end > $max)) {
                $tick++;
            }
            elsif ($end < $max) {
                $tick++;
            }
        }
        if ($tick > 0) {
            die("ERROR: Sequence from ctg $ctg was used more than once, e.g. in region $start-$end\n");
        }
        $max = $end unless ($max && ($max > $end));
    }
}


##############################

exit(0) unless ($csvin);

# creating csv file and renamed fasta file

my %chrom;
my %ind;
my %done;
open(CSV,$csvin);
open(CSVOUT,">${newout}.inter.csv");
open(CSVOUU,">${newout}.chromosome.list.csv");
my $secondout = Bio::SeqIO->new(-format => 'fasta',
			     -file => "> ${newout}.curated_primary.no_mt.unscrubbed.fa");


my $loc;
#open(TSVC,">${out}.tsv");
#open(TSVU,">${out}.unloc.tsv");

# gathering intel

while (<CSV>) {
    chomp;
    $loc++;
    my $entryname = $loc;
    my @a = split ",";
    my $total;
    my %part;
    foreach my $a (@a) {
        #print STDERR "searching for $a with $created{$a}\n";
        $total += $created{$a};   
        $part{$a} = $created{$a};          
        $entryname = "Z" if ($a =~ /Z/);
        $entryname = "W" if ($a =~ /W/);
        $entryname = "X" if ($a =~ /X/);
        $entryname = "Y" if ($a =~ /Y/);
    } 
    $chrom{$entryname}->{total} = $total;
    @{$chrom{$entryname}->{parts}} = (reverse sort {$part{$a} <=> $part{$b}} keys %part);
    
}   

# sorting and naming by size

my $number;
my $seenx;
foreach my $c (reverse sort {$chrom{$a}->{total} <=> $chrom{$b}->{total}} keys %chrom) {
    my $chrname;
    #print STDERR "checking $c\n";
    if ($c =~ /Z|W|Y|X/) {
        $chrname = $c;
        #print STDERR "found $c\n";
        
        # check that X/Z are bigger than Y/W
        if ($c =~ /Y|W/) {
            warn("Error: Sex chromosomes might be the wrong way around as found Y/W not smaller than X/Z\n") unless ($seenx);
        }
        $seenx++ if ($c =~ /Z|X/);
    }
    
    else {
        $number++;
        $chrname = $number;
    }
    my $first = 1;
    my $no = 1;
    foreach my $sub (@{$chrom{$c}->{parts}}) {

	die "ERROR: $sub (from $csvin) is not in the TPF ($tpf)\n " unless $assembly{$sub};

        if ($first) {
            print CSVOUT "$sub,$chrname,yes\n";
#            print TSVC   "$sub\t$chrname\tChromosome\n";
            print CSVOUU "SUPER_${chrname},$chrname,yes\n";
            my $newobj = $assembly{$sub};
            $newobj->display_id("SUPER_${chrname}");
            $secondout->write_seq($newobj);
            $done{$sub}++;
        }
        else {
            print CSVOUT "$sub,$chrname,no\n";
#            print TSVU   "$sub\t$chrname\n";
            print CSVOUU "SUPER_${chrname}_unloc_$no,$chrname,no\n";
            my $newobj = $assembly{$sub};
            $newobj->display_id("SUPER_${chrname}_unloc_$no");
            $secondout->write_seq($newobj);
            $done{$sub}++;
            $no++;
        }
        undef $first;
    }
}

foreach my $key (keys %assembly) {
    $secondout->write_seq($assembly{$key}) unless (exists $done{$key});
}

