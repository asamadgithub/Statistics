#!/usr/bin/perl

# reads xyz file and make all coordinates periodic
# yang@mpie.de
# 05.18.2016

print "put in \$alat first \n";
$alat[1]=54.06949667;
$alat[2]=26.75746154;
$alat[3]=25;
#print "$alat[1]\n";

open(IN,$ARGV[0]);
@in=<IN>;
$Nlines=$.;
close(IN);
#print "$Nlines\n";

($Nmol) = split(" ",$in[0]); # number of atoms
$Nrun=$Nlines/($Nmol+2);   # number of stemps
print "$Nmol, $Nrun\n";

# read data
for ($step=0;$step<$Nrun;$step++) {
  for ($mol=0;$mol<$Nmol;$mol++) {
    $line = ($Nmol+2)*$step+$mol+2;
    ($new[$step][$mol][0],$new[$step][$mol][1],$new[$step][$mol][2],$new[$step][$mol][3]) = split(" ",$in[$line]);
    # processing
    for ($i=0;$i<10;$i++) { 
       for ($k=1;$k<4;$k++) {
           if ($new[$step][$mol][$k] < 0) {
           #print "$new[$step][$mol][$k] \n";
           $new[$step][$mol][$k]=$new[$step][$mol][$k]+$alat[$k];            
           } elsif ($new[$step][$mol][$k] > $alat[$k]) {
           $new[$step][$mol][$k]=$new[$step][$mol][$k]-$alat[$k]; 
           }
       }
    }
  }
}



# write output data in new xyz
open(OUT,">$ARGV[0]2.xyz");
for ($step=0;$step<$Nrun;$step++) {
  print OUT $in[0];
  print OUT "Lattice=\"$alat[1] 0.0 0.0 0.0 $alat[2] 0.0 0.0 0.0 $alat[3]\" \n";
  for ($mol=0;$mol<$Nmol;$mol++) {
    print OUT "$new[$step][$mol][0] $new[$step][$mol][1] $new[$step][$mol][2] $new[$step][$mol][3] \n";
  }
#  print OUT "\n";
}
close(OUT);







