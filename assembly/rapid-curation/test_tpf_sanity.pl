#!/software/bin/perl -w

=head1 Name: test_tpf_sanity.pl

A script to run some sanity test on a TPF file

=head1 Usage: perl tpf_sanity.pl [-scafflevel] FILE

FILE is a file in TPF format. Either contig or scaffold level.

=head2 Options:

=over 12

=item -scafflevel

to use on scaffold level TPFs. Turns onchecks for naming, order and gaps 

=back

=cut

use strict;
use Getopt::Long;

my $scafflevel;
GetOptions(
    'scafflevel' => \$scafflevel, # if running on scaffold level TPF
)||die(`perldoc $0`);

# READ TPF
my $tpf = shift;
$tpf||die("missing FILE\n",`perldoc $0`);
open(TPF,$tpf);
my (%ctg,%scaff,%seen);

my $afterscaffgap;
my $afterctggap;
my $lastscaff;
my $lastline;
my %scaffseen;
my %fragCoords;

while (<TPF>) {
    my $line = $_;
    chomp($line);
    
    
    my @a=split;
    # non gaps
    if (/^\?/) {
        my ($ctg,$scaff) = ($a[1],$a[2]);
       
	print "line:$. please use an scaffold id in the form of xyz_123\n" unless $scaff = ~/_\d+$/;
 
        if ($lastscaff && ($lastscaff ne $scaff)) {
            if (exists $scaffseen{$scaff}) {
                print "$line\tThis scaffold name has been used before, please check and rename\n";
            }
        }
        $scaffseen{$scaff}++;
        


	# extract parent and coordinates to check Rapid Curation fragments for overlaps
	if ($ctg=~/([\w_]+):(\d+)-(\d+)/){
		if ($fragCoords{$1}) {
			my @matches;
			if (@matches = grep {$_->{start} <= $3 && $_->{end} >= $2} @{$fragCoords{$1}}){
				print "WARNING: overlapping coordinates for $1:$2-$3\n";
			        map {print "with $1:${\$_->{start}}-${\$_->{end}}\n"} @matches;	
			}
		}else{
			$fragCoords{$1}=();
		}
		push @{$fragCoords{$1}},{start => "$2",end => "$3"}
	}



        print "$line\tIs this an incorrectly retained haplotig?\n" if (($ctg =~ /_H/) || ($scaff =~ /_H/));
        print "$line\tIncorrect number of elements in line, did you forget the PLUS?\n" unless (scalar @a == 4);
        print "$line\torientation incorrect\n" unless (($a[3] eq "PLUS") || ($a[3] eq "MINUS"));
        # print "$line\tContig and scaffold have the same name\n" if ($ctg eq $scaff);
        
        # remember for later
        $ctg{$ctg}++;
        $scaff{$scaff}++;
        if ($lastscaff) {
            if ($afterscaffgap) {
                undef $afterscaffgap;
                print "$line\tnew scaffold name needed after scaffold gap\n" if ($lastscaff eq $scaff);
            }
            elsif ($afterctggap) {
                undef $afterctggap;
                print "$line\tsame scaffold name needed after contig gap\n" if ($lastscaff ne $scaff);
            }
            else {
                print "$line\tscaffold name incorrect\n" unless (($lastscaff eq $scaff) || $scafflevel);
		check_missing_gap($lastline,@a) if ($lastscaff eq $scaff);
            }
        }
        $lastscaff = $scaff;
        if ($lastline && $lastline !~/GAP/){
             check_orientation($lastline,@a);
        }
    }
    
    # gaps
    else {
        print "Empty line (this will generate two \"typo\" warnings)\n" unless (/\S+/);
        print "$line\ttypo\n" unless (/^GAP\s+TYPE-\d\s+\d+/i);
        $afterscaffgap++ if (/TYPE-3/i);
        $afterctggap++ if (/TYPE-2/i);
	print "duplicated line (at line $.): $_" if $lastline eq $_;
    }

    $lastline=$_;
}    

# Count occurrences
foreach my $ctg (keys %ctg) {
    if ($ctg =~ /(\S+)_\d+/) {
        if (exists $ctg{$1}) {
            print "$ctg\tunbroken retained: $1\n"; 
        }
    }
    if ($ctg{$ctg} > 1) {
        print "$ctg duplicated\n";
    }
}

# will be only run if there was no gap in front of it due to the call placement
# warn if it is the same scaffold, but different contigs without a gap in between
sub check_missing_gap{
	my $last = shift @_;
	my @currentLine = @_;
	my @lastLine = split /\s/, $last;

	my $l = $_; # copy the current line in
	chomp $l;

	if ($currentLine[2] eq $lastLine[2]){
	   my $currentId = "$1" if $currentLine[1]=~/(scaffold_\d+)/;
	   my $oldId = "$1"if  $lastLine[1]=~/(scaffold_\d+)/;

	   # for cases of split sequences
	   if ($currentLine[1] =~ /\d+\.\d+_\d+$/ && $lastLine[1] =~ /\d+\.\d+_\d+$/){
		    my $current = "$1" if $currentLine[1]=~ /(scaffold_\d+\.\d+_)\d+/;
	        my $old = "$1" if $lastLine[1]=~ /(scaffold_\d+\.\d+_)\d+/;
		    print "$l\t- GAP feature potentially missing\n" if ($current && $old) && ($current eq $old);
	   }

	   if ($currentLine[1]=~/scaffold_\d+\:\d+\-\d+/){ # as in: it is a rapid curation TPF
                	$currentId = "$1" if $currentLine[1]=~/(scaffold_\d+\:\d+\-\d+)/;
                	$oldId = "$1"if  $lastLine[1]=~/(scaffold_\d+\:\d+\-\d+)/;
           }
	   print "$l\t- GAP feature potentially missing\n" if ($currentId && $oldId )&& ($currentId ne $oldId);

	   # scaffold_1.14_1
	   # scaffold_1.14_3
	   # /scaffold_[XYZW\d]\.\d+_\d+/
   }  
	
}

# will only work for scaffold_\d+ style ids
sub check_orientation{
	my $last = shift @_;
	my @currentLine = @_;
	my @lastLine = split /\s/, $last;

        my $l = $_; # copy the current line in
	chomp $l;

	# skip the ones where the scaffold is different
	return unless ($lastLine[2] eq $currentLine[2]);

	# ignore anything where the regex doesn't match
	if ($currentLine[1]=~/(scaffold_\d+)\.(\d+)/){
	        my $id1 = "$1";
		my $count ="$2";
		if ($lastLine[1]=~/(scaffold_\d+)\.(\d+)/){
		    my $oldCount  = "$2";
	            my $id2 = "$1";
        	    return unless $id1 eq $id2;
		    if ($oldCount < $count && ($currentLine[3] eq 'MINUS' || $lastLine[3] eq 'MINUS')){ print "$l\t- MINUS set, but ascending IDs for $currentLine[1] and $lastLine[1]\n"}
		    if ($oldCount > $count && ($currentLine[3] eq 'PLUS' || $lastLine[3] eq 'PLUS')) { print "$l\t- PLUS set, but descending IDs for $currentLine[1] and $lastLine[1] ($count / $oldCount)\n"}
		}
	}
	

}
