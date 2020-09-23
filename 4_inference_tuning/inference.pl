#!/usr/bin/perl
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







#step_soliton_function();

phase_boundry();




sub phase_boundry{
    from_right_chiral_soliton();
} #end, chiral_soliton()



sub from_right_chiral_soliton{
    
    # 1-reading all relevent files
    #-----------------------------
    my ($ref_final_array,$n_atoms,$ref_hash_coeff,$ref_array_modes,$ref_file_array_1) = reading_files();
    my @final = @{$ref_final_array}; my %hash_ideals;  my @file_array_1 = @{$ref_file_array_1};
    # lattice constant and phonon mode application variables
    # -----------------------------------------------------
    my $stefan_a0= 3.86210691;my $b0= 25.64257;my $c0= 40;
    my $cp=13;my $z; #my $L= (2*$cp*(2*$stefan_a0))/4;
    my $l=$cp*2; my $ll=(4+2*($cp-1));
    my $L= (3+2*($cp-1))*$stefan_a0;               # defining period
    # ------------------------------------------------
    # creating new directory named "2__right_chiral_processing"
    # ------------------------------------------------
    my $dir_main=cwd;
    mkdir "$dir_main/gradual_PH_boundry",0755;
    chdir  "$dir_main/gradual_PH_boundry";
    my $dir_rc = cwd;
    
    #```````````````````````````````````````````````````````````````````````````````````````````
    #```````````````````````````````````````````````````````````````````````````````````````````
    # dynamical # of copies of ideal (8x2) structure in %hash_ideals
    # --------------------------------------------------------------
    foreach my $i(1..$cp){
        push @{ $hash_ideals{$i} }, @{dclone (\@final)};
    }
    
    # translating of ideal (8x2) structure in %hash_ideals
    # --------------------------------------------------------------
    foreach my $i(1..$cp){
        foreach my $atom(0..$#{$hash_ideals{$i}}){
            ${$hash_ideals{$i}}[$atom][0]+=(2*$stefan_a0*($i-1));
        }
    }

    # particular mode application
    # --------------------------------------------------------------
    my @right_rc_half; my @units_rc_half; my $ssslar_rc_half;my %z4_chiral_rc_half; $z=15;  ###############################
    foreach my $i(1..$cp){
        (undef,undef,$right_rc_half[$i],$ssslar_rc_half) = sending($i,\$hash_ideals{$i},$n_atoms,$ref_hash_coeff,$ref_array_modes,$ref_file_array_1,$z,$L);
        my @ff_rc_half =@{$right_rc_half[$i]};
        push @{$z4_chiral_rc_half{$i}}, @{dclone (\@ff_rc_half)};
        $units_rc_half[$i]=trimming_for_right_chiral(\@{$z4_chiral_rc_half{$i}},$i);
        
    }


    my $n3=$cp*104;
    open(my $fh_OUT_2 , "> PH_boundry_$l.xyz");
    print $fh_OUT_2 "$n3\n\n";
    foreach my $i(1..$cp){
        foreach my $atom(0..103){
            foreach my $k(0..3){
                print $fh_OUT_2 "$units_rc_half[$i][$atom][$k]\t";
            }print $fh_OUT_2 "\n";
        }
    }
    close ($fh_OUT_2);
    my $var_PH_bound = file_2d_array("PH_boundry_$l.xyz");
    third_subs_layer($l,$var_PH_bound,$dir_main,$dir_rc,$stefan_a0);
    `mv e__combine_as_2sub.xyz gradual_$l.xyz`;
    
    
    
    
    # two POSCAR file, old format with xyz2poscar.pl script
    # -----------------------------------------------------
    my $x_first=$stefan_a0*$l;
    `xyz2poscar.pl gradual_$l.xyz $x_first,0,0:0,$b0,0:0,0,$c0 > P_$l`;
    
    # writing gradual phase  in proper poscar format
    # -----------------------------------------------
    my $var_p2=file_2d_array("P_$l");
    my $file_2 = poscar_modifier($var_p2,$l); my @gradual=@{$file_2};
    
    
    # ppppp
    # -----
    
    open(my $fh_OUT3 , "> POSCAR_gradual_$l");
    foreach my $i(0..$#gradual){
        my $boxes = scalar @{$gradual[$i]};
        foreach my $j(0..($boxes-1)){
            print $fh_OUT3 "$gradual[$i][$j] ";
        }print $fh_OUT3 "\n";
    }
    close ($fh_OUT3);


    
} # end subroutine right_chiral_soliton();


sub trimming_for_right_chiral{
    my $var_rc = $_[0];  my  $unit_number = $_[1];
    my @un_trim_rc= @{$var_rc};
    my @rc_trim; chomp @rc_trim; my $trim_offset;
    $trim_offset = $un_trim_rc[158][1]+1;
    #----------------------------------------------------
    my $n1=104;
    #open (my $fh_OUT_1, "> local_$unit_number.xyz");
    #print $fh_OUT_1 "$n1\n\n";
    foreach my $i(0..$#un_trim_rc){
        if ($un_trim_rc[$i][1] < ($trim_offset)){
            foreach my $j(0..3){
                $rc_trim[$i][$j]=$un_trim_rc[$i][$j];
                #print $fh_OUT_1 "$rc_trim[$i][$j]\t";
            }#print $fh_OUT_1 "\n";
        }
    }
    #close ($fh_OUT_1);
    @rc_trim = grep { defined && m/[^\s]/ } @rc_trim;
    my $xx= scalar @rc_trim; #print "$xx\n";
    return (\@rc_trim);
    
    
} # end, sub trimming_for_right_chiral

#**************************************************************************************************************************
#       # modification in POSCAR
#**************************************************************************************************************************

sub poscar_modifier{
    my @pos = @{$_[0]};  my $size= $_[1];
    @pos = grep { defined && m/[^\s]/ } @pos;
    
    my (@tt1,@tt2); @tt2 = ("Si","In","H");
    foreach my $i(0..5){
        if ($i<3){
            $tt1[0][$i]=$pos[4][$i];
            #print "$tt1[0][$i] ";
        }#print "\n";
    }
    @tt1 = grep { defined && m/[^\s]/ } @tt1;
    
    my @title=("This is step soliton of size 8x$size");
    splice @pos,0,0,\@title;
    splice @pos,5,1;
    splice @pos,5,0,\@tt2;
    splice @pos,6,0,@tt1;
  
    return \@pos;
}

#**************************************************************************************************************************
    # Subrotutine to make structure with 3 substrate layers
#**************************************************************************************************************************
sub third_subs_layer{
    my $subs_size=$_[0];  my @prev_structure = @{$_[1]}; my $dir_main= $_[2];  my $dir_current = $_[3]; my $stefan_a0 = $_[4];
    my @geom = @{dclone(\@prev_structure)};
    copy("$dir_main/substrate.xyz","$dir_current/") or die "Copy failed: $!";
    my $unit_subs = file_2d_array('substrate.xyz');  my @array_unit_subs = @{$unit_subs};
    splice @array_unit_subs,0,2;
    
    
    # generateing substrate of required size
    # ---------------------------------------
    my %new_substrate;
    foreach my $i(1..$subs_size){
        push @{  $new_substrate{$i}  },  @{dclone(\@array_unit_subs)};
            foreach my $atom(0..$#{$new_substrate{$i}}){
                ${new_substrate{$i}}[$atom][1]+=($stefan_a0*($i-1));
            }
        }
    
    # ppppppppppppppppppppppppppppppppp
    # ---------------------------------
    my $n1=60*$subs_size;
    open(my $fh_OUT_1 , "> d__subs_3layer.xyz");
    print $fh_OUT_1 "$n1\n\n";
    foreach my $i(1..$subs_size){
        foreach my $atom(0..$#{$new_substrate{$i}}){
            foreach my $k(0..3){
                print $fh_OUT_1 "${$new_substrate{$i}}[$atom][$k]\t";
            }print $fh_OUT_1 "\n";
        }
    }
    close ($fh_OUT_1);
    my $var_total_substrate = file_2d_array('d__subs_3layer.xyz'); my @array_total_subs = @{$var_total_substrate};
    
    # extracting only In atoms from two laywered 'prev_structure'
    #------------------------------------------------------------
    my @in_only;
    foreach my $i(1..$#prev_structure){
        if ($prev_structure[$i][0] =~ m/^In/){
                foreach my $j(0..3){
                    $in_only[$i][$j]=$prev_structure[$i][$j];
                } $in_only[$i][3]-= 5;
            }
        }
    @in_only = grep { defined && m/[^\s]/ } @in_only;
    my $in_atoms= scalar @in_only;
    
    # extracting only Si atoms in two laywered soliton
    #-------------------------------------------------
    my @si_only;
    foreach my $i(1..$#array_total_subs){
        if ($array_total_subs[$i][0] =~ m/^Si/){
            foreach my $j(0..3){
                $si_only[$i][$j]=$array_total_subs[$i][$j];
            }
        }
    }
    @si_only = grep { defined && m/[^\s]/ } @si_only;
    my $si_atoms= scalar @si_only;
    
    # extracting only h atoms in two laywered soliton
    #-------------------------------------------------
    my @h_only;
    foreach my $i(1..$#array_total_subs){
        if ($array_total_subs[$i][0] =~ m/^H/){
            foreach my $j(0..3){
                $h_only[$i][$j]=$array_total_subs[$i][$j];
            }
        }
    }
    @h_only = grep { defined && m/[^\s]/ } @h_only;
    my $h_atoms= scalar @h_only;
    
    # combination as Si, In , H format
    # ----------------------------------------------
    my $total_atoms = $in_atoms+ $si_atoms+ $h_atoms;
    my @combined_modify;
    @combined_modify = (@si_only,@in_only,@h_only);
    
    # ppppppppppppppppppppppppppppppppp
    # ---------------------------------
    open(my $fh_OUT_2 , "> e__combine_as_2sub.xyz");
    print $fh_OUT_2 "$total_atoms\n\n";
    foreach my $i(0..$#combined_modify){
        foreach my $j(0..3){
            if($combined_modify[$i][1]<0){
                $combined_modify[$i][1]=0;
            }
             print $fh_OUT_2 "$combined_modify[$i][$j]\t";
        }print $fh_OUT_2 "\n";
    }
    close($fh_OUT_2);
    
    
    
    return (\@combined_modify);
} # end, subroutine third_subs_layer

#**************************************************************************************************************************
    # Subrotuine for read anyfile and return a 2d array
#**************************************************************************************************************************

 # end, sub step
 # -------------
sub file_2d_array{
    my $file= $_[0];
    my @file_array;
    open (my $fh_IN,"<",$file) or die "Could not open the $file.$!";
    while (my $line = <$fh_IN>){ chomp $line; $line=~ s/^\s+//g;
        if ($line =~ /^\w/){
            my @line_array=split (/\s+/,$line);
            push @file_array,\@line_array;
        }
    }
    
    return (\@file_array);
}  #end, sub reading_file


#**************************************************************************************************************************
    # Sub rotuine to apply phonon modes on a particular section of structure ,
#**************************************************************************************************************************
sub sending{
    my $key = $_[0]; my $recep = $_[1];my $n_atoms = $_[2]; my %coeff = %{$_[3]};my @modes = @{$_[4]};
    my @file_array_1 = @{$_[5]}; my $z= $_[6]; my $L=$_[7];
    my $pi = 3.14159265358979;
    # defining locally scoped variables inside sending
    # ------------------------------------------------
    my @final=@{$$recep};
    my @p;  @p =  @{dclone (\@final)};
    my @apply_mode=(520,526,517);         # 517 --> shear mode ; 520 & 526 --> rotary modes
    my @bottom_wire=(145,146,147,148,153,154,155,156,161,162,163,164,169,170,171,172);

    
    
    # logic to apply phonon mode
    # --------------------------
     chomp %coeff; chomp $z;
        if ($z<5){$coeff{517}=$coeff{517}} elsif(($z>5)&&($z<9)){$coeff{517}=-$coeff{517}}
    #if ($z==11){$coeff{517}=$coeff{517}}
        foreach my $i(@apply_mode){
            foreach my $atom(1..$n_atoms){
                foreach my $k(3..5){
                    next if (any {$atom eq  $_} @bottom_wire);
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                    if ($z eq 1){
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-1-> 520=NEG,526=NEG,
                    }
                    if ($z eq 2){
                        $coeff{526}=-$coeff{526};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-2-> 520=NEG,526=POS,
                        $coeff{526}=-$coeff{526};
                    }
                    if ($z eq 3){
                        $coeff{520}=-$coeff{520};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-3-> 520=POS,526=NEG,
                        $coeff{520}=-$coeff{520};
                    }
                    if ($z eq 4){
                        $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                        $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                    }
                    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                    if ($z eq 5){
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-5-> 520=NEG,526=NEG,
                    }
                    if ($z eq 6){
                        $coeff{526}=-$coeff{526};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-6-> 520=NEG,526=POS,
                        $coeff{526}=-$coeff{526};
                    }
                    if ($z eq 7){
                        $coeff{520}=-$coeff{520};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-7-> 520=POS,526=NEG,
                        $coeff{520}=-$coeff{520};
                    }
                    if ($z eq 8){
                        $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                        $p[$atom][$k-3]+=($coeff{$i}*$modes[$i][$atom][$k]);           # shear -ve;case-8-> 520=POS,526=POS,
                        $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
                    # ------------------------------------------------------------------------------------------------------------------
                    # 15+ ---> stands for the very gradual phase flip with pie periodic length
                    # ------------------------------------------------------------------------------------------------------------------
                    if ($z eq 15){
                        my $a= (cos(($p[$atom][0]/$L)*($pi)));   #print "$a\n";
                        
                        if ($a>0){
                            # This case is fully reduced to #case #4
                            if (($k==3) && ($p[$atom][1]>19)){
                                $coeff{517}=$coeff{517};
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                $p[$atom][$k-3]+=  $a * ($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};                 # modulation for uper chain, x-coordinates
                            }
                            elsif (($k!=3) && ($p[$atom][1]>19)){
                                $coeff{517}=$coeff{517};
                                $coeff{520}=-$coeff{520};  $coeff{526}=-$coeff{526};
                                $p[$atom][$k-3]+=  ($coeff{$i}*$modes[$i][$atom][$k]);
                                $coeff{520}=-$coeff{520};     $coeff{526}=-$coeff{526};
                            }

                            if (($k==3) && ($p[$atom][1]<19)){
                                
                                $coeff{517}=$coeff{517};
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                $p[$atom][$k-3]+=  $a* ($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};                 # modulation for uper chain, x-coordinates
                            }
                            
                            if (($k!=3) && ($p[$atom][1]<19)){
                                
                                $coeff{517}=$coeff{517};
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                $p[$atom][$k-3]+=   ($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};                 # modulation for uper chain, x-coordinates
                            }

=pod                            if ($p[$atom][1]<19){
                                
                                $coeff{517}=$coeff{517};
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};
                                $p[$atom][$k-3]+=   ($coeff{$i}*$modes[$i][$atom][$k]);           # shear +ve;case-4-> 520=POS,526=POS,
                                $coeff{520}=-$coeff{520};   $coeff{526}=-$coeff{526};                 # modulation for uper chain, x-coordinates
                            }
=cut
                            
                            
                            
                        }
                    }
                }
            }
        }

    
    # adjusting first column with relevent speice atoms for .xyz format
    # -----------------------------------------------------------------
    my @un_trim;chomp @un_trim;
    foreach my $atom(1..$n_atoms){
    ($un_trim[$atom][0],$un_trim[$atom][1],$un_trim[$atom][2],$un_trim[$atom][3])   =   ($file_array_1[$atom][0],$p[$atom][0],$p[$atom][1],$p[$atom][2]);
        #print  "$un_trim[$atom][0]---$un_trim[$atom][1]---$un_trim[$atom][2]---$un_trim[$atom][3]\n";
    }


  # selecting 8x2 from the processed 8x4 size
  # ------------------------------------------
    my @f_trim; chomp @f_trim;
    foreach my $i(1..$#un_trim){
        if ($un_trim[$i][1]<7.1){
            foreach my $j(0..3){
                $f_trim[$i][$j]=$un_trim[$i][$j];
                }
            }
        }

    @f_trim = grep { defined && m/[^\s]/ } @f_trim;
    @un_trim = grep { defined && m/[^\s]/ } @un_trim;
    my $ssslar1 = scalar @f_trim; my $ssslar2 = scalar @un_trim;
    
    
    @p=();
    return (\@f_trim,$ssslar1,\@un_trim,$ssslar2);    # @f_trim= trimmed structure, scalar @f_trim ,@un_trim)
}# end, sub sending

#**************************************************************************************************************************
    # Basic subroutine to read all input files and store the data in described arrays...
#**************************************************************************************************************************

sub reading_files{
    ## Useful Global Variables in this subroutine scope
    ##-------------------------------------------------
    my @all_files=("POSCAR","8x4_target.xyz","MODES","lincom.dat","POTCAR");      # 5 Files are stored in this array
    my @file_array_0;my $n_atoms;                                                 # To deal with "POSCAR" file,Total atoms in POSCAR
    my @file_array_1;                                                             # To deal with "8x4_target.xyz" geometry
    my @file_array_2; my @modes; my $n_modes=528;                                 # To deal with "MODES" file
    my @file_array_3 ;                                                        #  storing POTCAR file
    my %coeff;                                                                    # To deal with "lincom.dat" file
    
    
    my @equilibrium; my @final_1;
    my @evr;
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
                    $final_1[$atom][$k]= $modes[$mode_num][$atom][$k];
                    #print "$final[$atom][$k]\t";
                }#print "\n";
            }
        }
        $delta+=$n_atoms+2;
    }
    
    
    return (\@final_1,$n_atoms,\%coeff,\@modes,\@file_array_1);   #@file_array_1 , => target.xyz
    
    
} #end,  sub reading_file
