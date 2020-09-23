#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use Math::BigFloat;
use Math::Trig;
use List::MoreUtils qw(any);
use Storable qw (dclone);
use Cwd qw(cwd);
use File::Copy;
use File::Slurp;
use Scalar::Util qw(looks_like_number);
use POSIX;






my $dir = apply_modes(); 




########################################################################################################
# Apply 2 rotary and single shear mode to 8x4 , "stefan's structure"..
########################################################################################################

sub apply_modes{
## Useful Global Variables in this subroutine scope
##-------------------------------------------------
my @all_files=("POSCAR","8x4_target.xyz","MODES","lincom.dat","POTCAR");      # 5 Files are stored in this array
my @file_array_0;my $n_atoms;                                                 # To deal with "POSCAR" file,Total atoms in POSCAR
my @file_array_1;                                                             # To deal with "8x4_target.xyz" geometry
my @file_array_2; my @modes; my $n_modes=528;                                 # To deal with "MODES" file
my @file_array_3 ;                                                        #  storing POTCAR file
my %coeff;                                                                    # To deal with "lincom.dat" file
my @apply_mode=(520,526,517);         # 517 --> shear mode ; 520 & 526 --> rotary modes
#my @apply_mode=(517);
my @equilibrium; my @final; my @bottom_wire=(145,146,147,148,153,154,155,156,161,162,163,164,169,170,171,172);
my @evr;
my @P_AA;my @P_BB;my @P_AB;my @P_BA;my @N_AA;my @N_BB;my @N_AB;my @N_BA;
#**********************************************************************************************

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
                 if ($i eq $all_files[4]){
                     push @file_array_3,\@line_array;
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
#**************************************************************************************************

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
#***************************************************************************************************

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

#***************************************************************************************************

## Apply specific modes and getting structure after that
## -----------------------------------------------------
## [ A= -negative ; B= +positive ]  &&  [ P= +POSITIVE ; N= -NEGATIVE
@P_AA =  @{dclone (\@final)}; @P_BB =  @{dclone (\@final)};
@P_AB =  @{dclone (\@final)}; @P_BA =  @{dclone (\@final)};
@N_AA =  @{dclone (\@final)}; @N_BB =  @{dclone (\@final)};
@N_AB =  @{dclone (\@final)}; @N_BA =  @{dclone (\@final)};
    
    
foreach my $shear(1..2){
    if ($shear eq 1){$coeff{517}=$coeff{517}} else {$coeff{517}=-$coeff{517}}
        foreach my $z(1..8){
            foreach my $i(@apply_mode){
                foreach my $atom(1..$n_atoms){
                    foreach my $k(3..5){
                        next if (any {$atom eq  $_} @bottom_wire);
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            if (($shear eq 1) && ($z eq 1)){
                                    $P_AA[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-1-> 520=NEG,526=NEG,
                            }
                            if (($shear eq 1) && ($z eq 2)){
                                    $coeff{526}=-$coeff{526};
                                    $P_AB[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-2-> 520=NEG,526=POS,
                                    $coeff{526}=-$coeff{526};
                            }
                            if (($shear eq 1) && ($z eq 3)){
                                    $coeff{520}=-$coeff{520};
                                    $P_BA[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-3-> 520=POS,526=NEG,
                                    $coeff{520}=-$coeff{520};
                            }
                            if (($shear eq 1) && ($z eq 4)){
                                    $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                    $P_BB[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                                    $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                            }
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                            if (($shear eq 2) && ($z eq 5)){
                                    $N_AA[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-5-> 520=NEG,526=NEG,
                            }
                            if (($shear eq 2) && ($z eq 6)){
                                    $coeff{526}=-$coeff{526};
                                    $N_AB[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-6-> 520=NEG,526=POS,
                                    $coeff{526}=-$coeff{526};
                            }
                            if (($shear eq 2) && ($z eq 7)){
                                    $coeff{520}=-$coeff{520};
                                    $N_BA[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-7-> 520=POS,526=NEG,
                                    $coeff{520}=-$coeff{520};
                            }
                            if (($shear eq 2) && ($z eq 8)){
                                    $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                    $N_BB[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-8-> 520=POS,526=POS,
                                    $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                            }
                        #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                        }
                    }
                }
            }
}


#************************************************************************************
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
    #printf "\n  Initial average distance to target structure: %.6f A\n\n",$dist1;
    #printf "\n  Final average distance to target structure: %.6f A\n\n",$dist2;

    
## -------------------------------------------------------------------------------------------------
    ## change directory for generating files.
    ##----------------------------------------------
    my $dir = cwd;
    #print "\n$dir\n\n\n";
    mkdir "$dir/1__alles_possible_structures",0755;
    chdir "$dir/1__alles_possible_structures";
    ##----------------------------------------------
    
my $fh_OUT;
open($fh_OUT , "> 0__input.xyz");
print $fh_OUT "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $equilibrium[1][$atom][0]\t$equilibrium[1][$atom][1]\t$equilibrium[1][$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 1__P_AA__T.xyz");    #case-1:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is P_AA\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $P_AA[$atom][0]\t$P_AA[$atom][1]\t$P_AA[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 2__P_AB.xyz");    #case-2:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is P_AB\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $P_AB[$atom][0]\t$P_AB[$atom][1]\t$P_AB[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 3__P_BA.xyz");    #case-3:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is P_BA\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $P_BA[$atom][0]\t$P_BA[$atom][1]\t$P_BA[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 4__P_BB__T.xyz");    #case-4:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is P_BB\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $P_BB[$atom][0]\t$P_BB[$atom][1]\t$P_BB[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 5__N_AA.xyz");    #case-5:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is N_AA\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $N_AA[$atom][0]\t$N_AA[$atom][1]\t$N_AA[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 6__N_AB__T.xyz");    #case-6:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is N_AB\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $N_AB[$atom][0]\t$N_AB[$atom][1]\t$N_AB[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 7__N_BA__T.xyz");    #case-7:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is N_BA\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $N_BA[$atom][0]\t$N_BA[$atom][1]\t$N_BA[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 8__N_BB.xyz");    #case-8:
print $fh_OUT "$n_atoms\n";
print $fh_OUT "This is N_BB\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $N_BB[$atom][0]\t$N_BB[$atom][1]\t$N_BB[$atom][2]\n";
}
close ($fh_OUT);
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
open($fh_OUT , "> 9__end__pure_hexa.xyz");
print $fh_OUT "$n_atoms\n\n";
foreach my $atom(1..$n_atoms){
    print $fh_OUT "$file_array_1[$atom][0]  $file_array_1[$atom][1]\t$file_array_1[$atom][2]\t$file_array_1[$atom][3]\n";
}
close ($fh_OUT);
## ----------------------------------------------------------------------------------------------

    return ($dir);
} # end, sub apply modes()




