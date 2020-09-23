#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Math::BigFloat;
use Math::Trig;
use List::MoreUtils qw(any);
use Storable qw (dclone);

## 520,526 are trimmer rotary modes
## 517 is the useful shear mode

## Useful Global Variables
##------------------------
my @all_files=("POSCAR","8x4_target.xyz","MODES","lincom.dat");      # 3 Files are stored in this array
my @file_array_0;my $n_atoms;                                        # To deal with "POSCAR" file,Total atoms in POSCAR
my @file_array_1;                                                    # To deal with "8x4_target.xyz" geometry
my @file_array_2; my @modes; my $n_modes=528;                        # To deal with "MODES" file
my %coeff;                                                           # To deal with "lincom.dat" file
my @apply_mode=(517);
my @equilibrium; my @final; my @bottom_wire=(145,146,147,148,153,154,155,156,161,162,163,164,169,170,171,172);
my @evr;
my @shear_1;my @shear_2;
#***********************************************************************************************************************

## Reading array @all_files , one by one
## -------------------------------------
foreach my $i(@all_files){
    open (my $fh_IN,"<",$i) or die "Could not open the $i.$!";
    while (my $line = <$fh_IN>){ chomp $line; $line=~ s/^\s+//g;
        if ($line =~ /^\w/){
            my @line_array=split (/\s+/,$line);
            
            if ($i eq $all_files[0]){
                push @file_array_0,\@line_array;
            }
            
            if ($i eq $all_files[1]){
                push @file_array_1,\@line_array;
            }
            
            if ($i eq $all_files[2]){
                push @file_array_2,\@line_array;
            }
            
            if ($i eq $all_files[3]){
                $coeff{$line_array[0]}=$line_array[1];          # Keys are modes#, and values are their coefficients
            }
        }
    }
}


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
    }
    $delta+=$n_atoms+2;
}
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************
## Apply shear mode twice ----> hexa_forbidden structure
## -----------------------------------------------------
@shear_1 =  @{dclone (\@file_array_1)}; @shear_2 =  @{dclone (\@file_array_1)};

foreach my $z(1..2){
    foreach my $i(@apply_mode){
        foreach my $atom(1..$n_atoms){
            foreach my $k(3..5){
                #next if (any {$atom eq  $_} @bottom_wire);
                if ($z eq 1){
                $shear_1[$atom][$k-2]-=(    $coeff{$i}  *   ( 2 *   $modes[$i][$atom][$k])  );          # case-1-> 517=POS;
                }
                
                if ($z eq 2){
                    $coeff{517}=-$coeff{517};
                    $shear_2[$atom][$k-2]-=(    $coeff{$i}  *   ( 2 *   $modes[$i][$atom][$k])  );          # case-1-> 517=POS;
                    $coeff{517}=-$coeff{517};
                }
            }
        }
    }
}
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************
#**************************************************************************************************************************

open(my $fh_OUT_1 , "> 1.xyz");
print $fh_OUT_1 "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT_1 "$shear_1[$atom][0]  $shear_1[$atom][1]\t$shear_1[$atom][2]\t$shear_1[$atom][3]\n";
}


#~~~~~~~~~~

open(my $fh_OUT_2 , "> 2.xyz");
print $fh_OUT_2 "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT_2 "$shear_2[$atom][0]  $shear_2[$atom][1]\t$shear_2[$atom][2]\t$shear_2[$atom][3]\n";
}

