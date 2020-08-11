#!/usr/bin/perl

die "Usage: perl peak-caller.pl <Short WIG> <Log2(Short/Nuc) WIG> <BED OUT>\n" if(!$ARGV[2]);

BEGIN { push @INC, 'lib' }
use ngs;

($tMseq,$tMarr) = &ngs::readwigArray($ARGV[0]);
%Mseq  = %{$tMseq};
@Mseq_ar = @{$tMarr};

$tl2 = &ngs::readwig($ARGV[1]);
%l2  = %{$tl2};


($chr_me,$chr_std,$number) = &mean_se(\@Mseq_ar);
$str=sprintf"ALL: %.3f\t%.3f\t%d\n",$chr_me,$chr_std,$number;
print $str;

@Mseq_ar = sort {$a <=> $b} @Mseq_ar;
print "Start: $Mseq_ar[$#Mseq_ar]\n";

$Npeaks=0;
$prev_peaks=-1;
$iter=0;
$cutoff=4;

while($Npeaks>$prev_peaks){
	print "$iter Npeak: $Npeaks Prev_peak: $prev_peaks N:$#Mseq_ar\n";
	$prev_peaks=$Npeaks;
	for($i=$#Mseq_ar;$i>=0;$i--) {
		$value = ($Mseq_ar[$i]-$chr_me)/$chr_std;
		if( $value >= $cutoff){
			splice(@Mseq_ar,$i,1);	
			$Npeaks++;
		}else{
			break;
		}
	}
	($chr_me,$chr_std,$number) = &mean_se(\@Mseq_ar);
	$iter++;
}
print "$iter Npeak: $Npeaks Prev_peak: $prev_peaks N:$#Mseq_ar\n";
print "\nEnd: $Mseq_ar[$#Mseq_ar] Mean: $chr_me SD: $chr_std\n\n";

$cutoff=$Mseq_ar[$#Mseq_ar];

print "Cutoff: $cutoff\n";
#($chr_me,$chr_std,$number) = &mean_se(\@Mseq_ar);
#$str=sprintf"ALL: %.3f\t%.3f\t%d\n",$chr_me,$chr_std,$number;
#print $str;

#die "test ends here\n";

open(BED,">$ARGV[2]") || die "$!\n";
foreach $i ( keys(%Mseq) ){ #each chromosome
	$#sorted_keys=-1;
	$#peak_array=-1;
	$#val_array=-1;
	@sorted_keys = sort {$a <=> $b} (keys(%{$Mseq{$i}}) );
	foreach $j (@sorted_keys){ # Find peaks
		$value = ($Mseq{$i}{$j}-$chr_me)/$chr_std ;
		if( $Mseq{$i}{$j} > $cutoff && $l2{$i}{$j} >= 2){ # Higher than iteratively estimated cut-off and higher than log2 of 2
			push(@peak_array,$j); 
			push(@val_array,$value);
		}
	}
	print STDERR "$i $#peak_array $#val_array\n";
	$#collapse=-1;	$#vcollapse=-1;	$prev=-1;
	for($k=0;$k<=$#peak_array;$k++){ # Find summits
		if($k==0 ){
			push(@collapse,$peak_array[$k]);
			push(@vcollapse,$val_array[$k]);
		}elsif($peak_array[$k]==$prev+1){
			push(@collapse,$peak_array[$k]);
			push(@vcollapse,$val_array[$k]);
		}elsif($#collapse>-1){
			($vpos,$vval)=&max_pos(\@vcollapse);
			#print "vpos: $vpos @vcollapse \n";
			$me=$collapse[$vpos];
			$me = $me;
			$number=$#collapse+1;
			print "Chr: $i Pos: $me Z-score: $vval numPoints: $number vpos: $vpos l2val: $l2{$i}{$me}\n"; #if($vval>=2 && $lme>=2); # peak-chromosome, peak-position, z-score at peak-position, number of points making up the peak, log2 value at peakif($
			if($number<2){
				$strt=$me-1;
				$end=$me+1;
			}else{
				$strt=$collapse[0];
				$end=$collapse[$#collapse];
			}
			$str=sprintf("%s\t%d\t%d\t%s_%d\t%5.1f\n",$i,$strt,$end,$i,$me,$vval); 
			print BED $str; 
			$#collapse=-1;
			$#vcollapse=-1;
			$#lcollapse=-1;
			push(@collapse,$peak_array[$k]);
			push(@vcollapse,$val_array[$k]);
		}
		$prev = $peak_array[$k];
	}
}
close(BED);

sub mean_se {
  my @ar = @{$_[0]};
  my $sum=0;
  foreach my $ii (@ar){
    $sum+=$ii;
  }
  my $mean=$sum/($#ar+1);
  $sum=0;
  foreach my $ii (@ar){
    $x=($ii-$mean)**2;
    $sum+=$x;
  }
  my $se=sqrt($sum/($#ar+1));
	my $num = $#ar+1;
  return ($mean,$se,$num);
}

sub max_pos {
	my @ar = @{$_[0]};
	my @sar = sort {$a<=>$b} @ar;
	my $pos;
	for(my $ii=0;$ii<=$#ar;$ii++){
	 if($ar[$ii] == $sar[$#sar]){
	 	$pos=$ii;
		last;
	 }
	}
	return ($pos,$sar[$#sar]);
}
