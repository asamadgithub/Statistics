#!/usr/bin/perl

=pod
Searches for the linear combination of eigenvectors in
MODES that minimizes the distance to the target structure
A.Samad, 20.06.2017
=cut

use strict;
use warnings;
use Data::Dumper;
use Math::BigFloat;
use Math::Trig;


## Useful Global Variables
##------------------------
my @all_files=("POSCAR","8x4_target.xyz","MODES","lincom.dat");      # 3 Files are stored in this array
my @file_array_0;my $n_atoms;                                        # To deal with "POSCAR" file,Total atoms in POSCAR
my @file_array_1;                                                    # To deal with "8x4_target.xyz" geometry
my @file_array_2; my @modes; my $n_modes=528;                        # To deal with "MODES" file
my %coeff;                                                           # To deal with "lincom.dat" file
my @apply_mode=(520,526);
my @equilibrium; my @final;
my @evr;
#***********************************************************************************************************************

## Reading array @all_files , one by one
## -------------------------------------
foreach my $i(@all_files){
  open (my $fh_IN,"<",$i) or die "Could not open the $i.$!";
     while (my $line = <$fh_IN>){ chomp $line; $line=~ s/^\s+//g;
         if ($line =~ /^\w/){
             my @line_array=split (/\s+/,$line);
             
                 if ($i eq $all_files[0]){ push @file_array_0,\@line_array;}
             
                 if ($i eq $all_files[1]){ push @file_array_1,\@line_array;}
             
                 if ($i eq $all_files[2]){push @file_array_2,\@line_array;}
             
                 if ($i eq $all_files[3]){
                     $coeff{$line_array[0]}=$line_array[1];          # Keys are modes#, and values are their coefficients
                 }
        }
    }
}
#print Dumper(\%coeff); #print "$coeff{526}\n";

                ##Total Number of atoms
                ##---------------------
                foreach my $j(0..$#{$file_array_0[6]}){
                        $n_atoms+=$file_array_0[6][$j];
                }#print "\n$n_atoms\n";
#*************************************************************************************************************************

## Module for fetching Geometries and Eigenvectors
## -----------------------------------------------
my $s_line=3; my $delta=0;

foreach my $mode_num(1..$n_modes){
        foreach my $j(($s_line+$delta)..($n_atoms+$s_line-1+$delta))  { my $atom=$j-$s_line-$delta+1;
            foreach my $k(0..5)      {
                $modes[$mode_num][$atom][$k]=$file_array_2[$j][$k];
        	}
            if($mode_num eq 1){
                foreach my $k(0..2)      {
                    $equilibrium[$mode_num][$atom][$k]=$modes[$mode_num][$atom][$k];
                    #print "$equilibrium[$mode_num][$atom][$k]\t";
                }#print "\n";
            }
            if($mode_num eq $n_modes){
                foreach my $k(0..2)      {
                    $final[$atom][$k]= $modes[$mode_num][$atom][$k];
                    #print "$final[$atom][$k]\t";
                }#print "\n";
            }
	 }
    $delta+=$n_atoms+2;
}
#**************************************************************************************************************************

=pod
## To readily looking for a particular mode
##  ---------------------------------------
 print "=========================================================================================\n";
foreach my $mode_num(50){    print "\nThis is mode number $mode_num f\n";
    print "\t            X\t            Y\t            Z\t            dx\t            dy\t            dz\n";
    foreach my $atom(1..$n_atoms){  print "atoms# $atom\t";
        foreach my $k(0..5){
            print "$modes[$mode_num][$atom][$k]\t";
            #print "[$mode_num][$atom][$k]\t";
        }print "\n";
    }
}
print "=========================================================================================\n";
=cut

#**************************************************************************************************************************

## Apply specific modes and getting structure after that
## -----------------------------------------------------
foreach my $i(@apply_mode){
    foreach my $atom(1..$n_atoms){
        foreach my $k(3..5){
#         print "$modes[$i][$atom][$k]\t"; print "$coeff{$i}"; print "$modes[$n_modes][$atom][$k-3]\t";
          $final[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);
#         print "$final[$atom][$k-3]\t";
        }#print "\n";
    }
}
#**************************************************************************************************************************

## -----------------------------------
## Section to write  & Display result
## -----------------------------------

##  Avg distances between atoms
##  --------------------------
my $dist1;my $dist2;
foreach my $atom(1..$n_atoms){
    foreach my $k(0..2){
        $dist1+=sqrt(($file_array_1[$atom][$k+1]-$equilibrium[1][$atom][$k])**2);
        $dist2+=sqrt(($file_array_1[$atom][$k+1]-$final[$atom][$k])**2);
    }
}
$dist1/=$n_atoms; $dist2/=$n_atoms;
printf "\n  Initial average distance to target structure: %.6f A\n\n",$dist1;
printf "\n  Final average distance to target structure: %.6f A\n\n",$dist2;

## ---------------------------------------------------------------------------------

open(my $fh_OUT , "> sam.xyz");

print $fh_OUT "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $equilibrium[1][$atom][0]\t$equilibrium[1][$atom][1]\t$equilibrium[1][$atom][2]\n";
}
print $fh_OUT "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $final[$atom][0]\t$final[$atom][1]\t$final[$atom][2]\n";
}

print $fh_OUT "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $file_array_1[$atom][1]\t$file_array_1[$atom][2]\t$file_array_1[$atom][3]\n";
}
close ($fh_OUT);
