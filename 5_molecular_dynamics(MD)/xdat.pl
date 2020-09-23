#!/usr/bin/perl
# Abdus-Samad, 16.01.18
# a.samad@mpie.de

print <<OUT;

Take XDARCAR file and calculate avg over all iterations
-------------------------------------------------------
Kindly remove the lines (0-6) for every subsequent XDATCAR file after first;
Only the first XDATCAR will contain such information;
e.g. cat XDATCAR1(nur lattice info) XDATCAR2 XDATCAR3 XDATCAR4 > XDATCAR

OUT

#------------------------------------------------------------
use strict;
use warnings;
use Data::Dumper;
use Math::BigFloat;
use Math::Trig;
use List::MoreUtils qw(any);  use List::Util 'reduce';
use Storable qw (dclone);
use Cwd qw(cwd);
use File::Copy;  use File::Slurp;
use Scalar::Util qw(looks_like_number);
use POSIX;
use Tie::IxHash;
use List::Util 'max';
#------------------------------------------------------------




my @all_files=('XDATCAR');
my @xyz_struc;
foreach my $i(@all_files){
    open (my $fh_IN,"<",$i) or die "Could not open the $i.$!";
    while (my $line = <$fh_IN>){ chomp $line; $line=~ s/^\s+//g;
        #if ($line =~ /^\w/){
            my @line_array=split (/\s+/,$line);
                if ($i eq $all_files[0]){
                    push @xyz_struc,\@line_array;
                }
        #}
    }
}
my $xyz_last_line= scalar @xyz_struc;
print "Total number of lines in XDATCAR file:                  $xyz_last_line\n";
my $var1= $xyz_struc[6][0]+$xyz_struc[6][1]+$xyz_struc[6][2];
my $last_itr_info=$xyz_last_line-$var1;
print "Line# containing info of last ieraation XDATCAR file:   $last_itr_info\n";


# Useful information about @AA extracted from @struc_xyz
# ------------------------------------------------------
my @AA; @AA =  @{dclone (\@xyz_struc)};           # Array @AA stats from first "Direct configuration=     1"
my $n_atoms= $AA[6][0]+$AA[6][1]+$AA[6][2];
my $natom_plus_1 = $n_atoms+1;                    # its # of atoms + 1
#print "$n_atoms\n";
#print "$natom_plus_1\n";
splice @AA,0,7;
my $lines= scalar @AA;                            # total number of lines exclusing first 6
$lines=$lines-$natom_plus_1;                      # it gives the line number of last_iteration wrt to  @AA
my $last_iteration = $AA[$lines][2];



#------------------------------------------------------------------------------
# Iteration according to desired operation over all iterations
#------------------------------------------------------------------------------
my @collect;
# its starting from 1 till no_ionic steps
#foreach my $iteration(1..$last_iteration){

foreach my $atom_num(1..$n_atoms){   #print "$atom_num:-->\n";
    foreach my $iteration(1..$last_iteration){
        my $v =(($iteration-1)*$natom_plus_1);
            foreach my $k(0..2){
                $collect[$atom_num][$k]+=$AA[$v + $atom_num][$k];
            }
    }
}
#print "@{$AA[409]}\n";                          # important checking line
#------------------------------------------------------------------------------
# Averaging the whole 2-d array
#------------------------------------------------------------------------------
my @FF; @FF =  @{dclone (\@collect)};
foreach my $i(1..$#FF){
    foreach my $j(0..$#{$FF[$i]}){
        $FF[$i][$j]/=$last_iteration;
        
    }
}

my $f= scalar @FF;

print "$last_iteration\n";
#------------------------------------------------------------------------------
# Printing of file
#------------------------------------------------------------------------------
open (my $fh_OUT, ">POSCAR_Avg");
foreach my $i(0..7){
    foreach my $j(0..$#{$xyz_struc[$i]}){
        print $fh_OUT "$xyz_struc[$i][$j] ";
    } print $fh_OUT "\n";}

foreach my $i(1..$#FF){
    foreach my $j(0..$#{$FF[$i]}){
        print $fh_OUT "$FF[$i][$j]    ";
    } print $fh_OUT "\n";}
close($fh_OUT);
print "POSCAR_Avg had bee written already\n\n";


