#!/usr/bin/perl
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

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Useful global variable
# ----------------------
my $nat; my $cfgs;

# Reading .xyz file for averaging  and POSCAR for file parameter
# --------------------------------------------------------------

my @all_files=('a.xyz');
my @xyz_struc;
foreach my $i(@all_files){
    open (my $fh_IN,"<",$i) or die "Could not open the $i.$!";
    while (my $line = <$fh_IN>){ chomp $line; $line=~ s/^\s+//g;
        if ($line =~ /^\w/){
            my @line_array=split (/\s+/,$line);
                if ($i eq $all_files[0]){
                    push @xyz_struc,\@line_array;
                }
            }
        }
}

$nat= $xyz_struc[0][0];
$cfgs= ((scalar @xyz_struc)/($nat+2))-1;    # -1 because the iteration numbers starts from zero
print "\n\nTotal # of atoms $nat\n";
my $cc=$cfgs+1;
print "Total # of ionic configuration $cc\n********************************************\n";

# Sum up all relevant co-ordinates
# --------------------------------
my $escape=1;                              # It will start from  $escape-1; <<<<<<<>>>>>>>>>>>>>>><<<<<<<<<<<<>>>>>>>>>>>
#                                          e.g. if $escape=2002, it will stating summing up from 2001 and onwards...
my @sum; my @avg;                           # summing up all relevant iterations
my $div=(($cfgs+1)-$escape)+1;              # this will use for divison after usmming up all

foreach my $atom_num(0..($nat-1)){
    foreach my $iteration($escape..($cfgs+1)){
        my $v =2+(($iteration-1)*($nat+2));
        print "$v  $atom_num\n";
        foreach my $k(1..3){
            $sum[$atom_num][$k]+=$xyz_struc[$v + $atom_num][$k];
        }
    }
}


# Slicing first column of @xyz_struc for atomic species
# -----------------------------------------------------
my @specie;
foreach my $i(2..($nat+1)){
    my $j=$i-2;
    $specie[$j][0]=$xyz_struc[$i][0];
    #print "$specie[$j][0]\n";
}

# Avergaing the whole data
# ------------------------
print "Divided by/Average of :$div ionic iterations\n";
foreach my $atom_num(0..($nat-1)){
    foreach my $k(1..3){
        $avg[$atom_num][$k]=($sum[$atom_num][$k]/$div);
        #print "$avg[$atom_num][$k]\t";
    }
    #print "\n";
}


# Writing POSCAR.xyz & POSCAR__avg
# --------------------------------
open (my $fh_OUT, ">avg.xyz");
print $fh_OUT "$nat\n";
print $fh_OUT "8x14 Avg Structure of last $div iterations\n";
foreach my $atom_num(0..($nat-1)){   print $fh_OUT "$specie[$atom_num][0]\t";
    foreach my $k(1..3){
        print $fh_OUT "$avg[$atom_num][$k]\t";
    }
    print $fh_OUT "\n";
}

`xyz2poscar.pl avg.xyz 54.06949667,0,0:0,26.75746154,0:0,0,25 > POSCAR___avg`;
print "avg.xyz and POSCAR__avg files has been written\n\n";



