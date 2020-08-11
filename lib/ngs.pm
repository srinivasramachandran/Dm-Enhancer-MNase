package ngs;
use strict;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);
$VERSION     = 1.00;

my %dm3_chrom_size = (
	"Uextra"	=> 29004656,
	"3R"			=> 27905053,
	"3L"			=> 24543557,
	"2L"			=> 23011544,
	"X"				=> 22422827,
	"2R"			=> 21146708,
	"U"				=> 10049037,
	"2RHet"		=>	3288761,
	"3LHet"		=>	2555491,
	"3RHet"		=>	2517507,
	"4"				=>	1351857,
	"2LHet"		=>	 368872,
	"YHet"		=>	 347038,
	"XHet"		=>	 204112,
	"M"				=>	  19517,
);

sub writeWig{
	#$jnk=writeWig(wig_hash,out_file,wig_step)
	my $href = $_[0];
	my $outfile = $_[1];
	my $step = 10;
	$step = $_[2] if($_[2] > 0);
	my ($href,$outfile)=@_;
	my %Mseq  = %{$href};
	open(OUT,">$outfile") || die "OUT writeWig $outfile $!\n";
	print OUT "track type=wiggle_0\n";
	foreach my $i (keys (%Mseq) ){
		print OUT "variableStep  chrom=chr$i span=$step\n";
		my %thash = %{$Mseq{$i}};
		foreach my $j ( sort {$a<=>$b} keys(%thash) ){
			print OUT "$j $Mseq{$i}{$j}\n" if($_[3] eq "dm3" && $j <= $dm3_chrom_size{$i});
		}
	}
	close(OUT);
	return(1);
}

sub writeWig_general{
	#$jnk=writeWig(wig_hash,out_file,wig_step)
	my $href = $_[0];
	my $outfile = $_[1];
	my $step = 10;
	$step = $_[2] if($_[2] > 0);
	my ($href,$outfile)=@_;
	my %Mseq  = %{$href};
	open(OUT,">$outfile") || die "OUT writeWig $outfile $!\n";
	print OUT "track type=wiggle_0\n";
	foreach my $i (keys (%Mseq) ){
		print OUT "variableStep  chrom=chr$i span=$step\n";
		my %thash = %{$Mseq{$i}};
		foreach my $j ( sort {$a<=>$b} keys(%thash) ){
			print OUT "$j $Mseq{$i}{$j}\n";
		}
	}
	close(OUT);
	return(1);
}


sub readYwig{
#wig format: Chr Coord Value
	my @ytemp; my %Ywig; my $ywig_line;	
	open(YWIG,$_[0]) || die "YWIG $!\n";
	while($ywig_line=<YWIG>){
		@ytemp = split /[\ \s\n\t]+/, $ywig_line;
		$ytemp[0]=~s/chr//;
		$Ywig{$ytemp[0]}{$ytemp[1]}=$ytemp[2]; # $Ywig{Chromosome}{coordinate} = value
	}
	close(YWIG);
	return (\%Ywig);
}

sub smoothwig_sub{
# outwig = smoothwig_sub(peak_file,wig_file,smooth_half_window,wig_step,peak_window)
	my @temp ;my %swig; my %wig; my $wig_line; my $chrid; my $step; my %pos;
	my %sm_flag;
	my $window = $_[2];
	my $step = $_[3];
	my $pwindow = $_[4];
	print "Smooth_window $window\n";
	print "Peak_window $pwindow\n";
	open(FILE,$_[0]) || die "Peak file $!\n";
	while(chomp(my $peak_line=<FILE>)){
		@temp = split /[\ \s\n\t]+/, $peak_line;
		$chrid = $temp[1];
		$chrid=~s/chr//;
		my $peak = (int($temp[3]/10+0.5))*10+1;
		for(my $i=$peak-$pwindow;$i<=$peak+$pwindow;$i+=$step){
			$pos{$chrid}{$i} = 1;
		}
	}
	close(FILE);
	my @kpos = keys(%pos);
	my $kpos_val = join("\t",@kpos);
	print "$kpos_val\n";
	open(WIG,$_[1]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
			print "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $wig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	
	foreach $chrid ( keys(%pos) ){
		print STDERR "New $chrid\n";
		my @sjunk = keys(%{$pos{$chrid}});
		print STDERR "$#sjunk\n";
		my @junk = sort {$a<=>$b} @sjunk ;
		print STDERR "$junk[$#junk]\n";
		my $i;
		my $j;
		foreach $i (@junk){
			my $sum=0;
			for($j=$i-$window;$j<=$i+$window;$j+=$step){
				$sum += $wig{$chrid}{$j};
			}
			$swig{$chrid}{$i}= $sum/($window*2+1) if( $sm_flag{$chrid}{$i} != 1 && $sum>0);
			$sm_flag{$chrid}{$i} = 1;
		}
	}
	return (\%swig);
}


sub smoothwig{
# outwig = smoothwig(wig_file,half_window,wig_step)
	my @temp ;my %swig; my %wig; my $wig_line; my $chrid; my $step;
	
	my $window = $_[1];
	
	my $step = $_[2];

	print "$window\n";
	
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
			print "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $wig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	
	foreach $chrid ( keys(%wig) ){
		print STDERR "New $chrid\n";
		my @sjunk = keys(%{$wig{$chrid}});
		print STDERR "$#sjunk\n";
		my @junk = sort {$a<=>$b} @sjunk ;
		print STDERR "$junk[$#junk]\n";
		my $sum=0;
		for(my $i = $junk[0]-$window;$i<= $junk[0]+$window;$i+=$step){
			if( exists $wig{$chrid}{$i}){
				$sum += $wig{$chrid}{$i};
			}
			$swig{$chrid}{$junk[0]}= $sum/($window*2+1);
		}
		for(my $j=$junk[0]+$step; $j <=$junk[$#junk];$j+=$step){
			$sum = $sum - $wig{$chrid}{$j-$window};
			$sum = $sum + $wig{$chrid}{$j+$window};
			$swig{$chrid}{$j}= $sum/($window*2+1);
		}
	}
	return (\%swig);
}


sub subwig{
# outwig = subwig(peakfile, window, wig_file, step)
	my @temp; my %wig; my $wig_line; my $chrid; my %pos;
	my $window = $_[1];
	my $step = 10;
	$step = $_[3] if($_[3] > 0);
	print "$window\n";
	open(FILE,$_[0]) || die "Peak file $!\n";
	while(chomp(my $peak_line=<FILE>)){
		@temp = split /[\ \s\n\t]+/, $peak_line;
		$chrid = $temp[1];
		$chrid=~s/chr//;
		my $peak = (int($temp[3]/10+0.5))*10+1;
		print "$chrid $peak\n";
		for(my $i=$peak-$window;$i<=$peak+$window;$i+=$step){
			$pos{$chrid}{$i} = 1;
		}
	}
	close(FILE);

	open(WIG,$_[2]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
			print "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1] if(exists $pos{$chrid}{$temp[0]}); # $wig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	return (\%wig);
}

sub make_peakhash {
	# hash_ref = make_peakhash(peak_file, window)
	my $chrid; my %pos; my $window = $_[1]; my $peak_file = $_[0];
	my @temp;
	open(FILE,$peak_file) || die "Peak file $!\n";
	while(chomp(my $peak_line=<FILE>)){
		@temp = split /[\ \s\n\t]+/, $peak_line;
		$chrid = $temp[1];
		$chrid=~s/chr//;
		my $peak = int(($temp[3]/10) + 0.5)*10 + 1;
		print "$chrid $peak\n";
		for(my $i=$peak-$window;$i<=$peak+$window;$i+=10){
			$pos{$chrid}{$i} = 1;
		}
	}
	print STDERR "Peak Filler Filled\n";
	close(FILE);
	return(\%pos);

}

sub subpairs{
# jnk = subpairs(peak_hash, pairs_file, out_file)
	my $href=$_[0];
	my %pos  = %{$href};
	my $pair_file = $_[1];
	my $outfile=$_[2];
	my @temp;
	print STDERR "pairs: $pair_file out_pairs: $outfile\n";
	open(OUT,'>>',$outfile) || die "OUT writeWig $outfile $!\n";
	open(PAIRS,$pair_file) || die "WIG $!\n";
	while(my $bed_line=<PAIRS>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $st = int(($temp[2]/10) + 0.5)*10 + 1;
		my $en = int(($temp[3]/10) + 0.5)*10 + 1;
		$temp[0]=~s/chr//;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif(exists $pos{$temp[0]}{$st} || exists $pos{$temp[0]}{$en}){
			print OUT $bed_line;
		}
	}
	close(PAIRS);
	close(OUT);
	print STDERR "Done with pairs file\n";
	return (1);
}

sub subbed{
# jnk = subpairs(peak_hash, pairs_file, out_file)
	my $href=$_[0];
	my %pos  = %{$href};
	my $bed_file = $_[1];
	my $outfile=$_[2];
	my @temp;
	print STDERR "BED: $bed_file out_bed: $outfile\n";
	open(OUT,'>>',$outfile) || die "OUT writeWig $outfile $!\n";
	if($bed_file=~/gz$/){
		open(BED,"gunzip -c $bed_file |") || die "cannot open pipe to $bed_file";
	}else{
			open(BED,$bed_file) || die "BED $!\n";
	}
	while(my $bed_line=<BED>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $st = int(($temp[1]/10) + 0.5)*10 + 1;
		my $en = int(($temp[2]/10) + 0.5)*10 + 1;
		$temp[0]=~s/chr//;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif(exists $pos{$temp[0]}{$st} || exists $pos{$temp[0]}{$en}){
			print OUT $bed_line;
		}
	}
	close(BED);
	close(OUT);
	print STDERR "Done with BED file\n";
	return (1);
}


sub readDwig{
#wig format: Chr Coord Value
	my @ytemp; my %Dwig; my $dwig_line;	
	open(DWIG,$_[0]) || die "YWIG $!\n";
	while($dwig_line=<DWIG>){
		@ytemp = split /[\ \s\n\t]+/, $dwig_line;
		$ytemp[0]=~s/chr//;
		$Dwig{$ytemp[0]}{$ytemp[1]}=$ytemp[3] if($ytemp[3]!=0); # $Dwig{Chromosome}{coordinate} = value
	}
	close(DWIG);
	return (\%Dwig);
}

sub readbedgraph{
#format: chr start_coord end_coord value
	my @temp; my %bedgraph; my $bed_line; my $chrid;
	open(BEDGRAPH,$_[0]) || die "WIG $!\n";
	while(my $bed_line=<BEDGRAPH>){
		@temp = split /[=\ \s\n\t]+/, $bed_line;
		if($bed_line=~/track/){
		}else{
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $st=( int( ($temp[1]+1) / 10 ) )*10+1;
			my $en=( int(  $temp[2]    / 10 ) )*10+1;
			for(my $i=$st;$i<=$en;$i+=10){
				$bedgraph{$chrid}{$i}=$temp[3];
			}
		}
	}
	close(BEDGRAPH);
	return (\%bedgraph);
}

sub readbedgraph50{
#format: chr start_coord end_coord value
	my @temp; my %bedgraph; my $bed_line; my $chrid;
	open(BEDGRAPH,$_[0]) || die "WIG $!\n";
	while(my $bed_line=<BEDGRAPH>){
		@temp = split /[=\ \s\n\t]+/, $bed_line;
		if($bed_line=~/track/){
		}else{
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $st=( int( ($temp[1]+1) / 50 ) )*50+1;
			my $en=( int(  $temp[2]    / 50 ) )*50+1;
			for(my $i=$st;$i<=$en;$i+=50){
				$bedgraph{$chrid}{$i}=$temp[3];
			}
		}
	}
	close(BEDGRAPH);
	return (\%bedgraph);
}



sub readbedgraph_bp{
#format: chr start_coord end_coord value
	my @temp; my %bedgraph; my $bed_line; my $chrid;
	open(BEDGRAPH,$_[0]) || die "WIG $!\n";
	while(my $bed_line=<BEDGRAPH>){
		@temp = split /[=\ \s\n\t]+/, $bed_line;
		if($bed_line=~/track/){
		}else{
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $st=$temp[1];
			my $en=$temp[2];
			for(my $i=$st;$i<=$en;$i++){
				$bedgraph{$chrid}{$i}=$temp[3];
			}
		}
	}
	close(BEDGRAPH);
	return (\%bedgraph);
}

sub readbed_c{
	my ($bed_file,$min,$max)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	open(BED,$bed_file) || die "WIG $!\n";
	while(my $bed_line=<BED>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $len = $temp[2]-$temp[1];
		if($#temp < 2){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif($len>=$min && $len<=$max){
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $mp = int( ($temp[1]+1+$temp[2])/2 +0.5);
			$twig{$chrid}{$mp}++;
			$ct++;
		}
	}
	close(BED);
	print STDERR "Lines used: $ct\n";
	my $GN = 139712364; # Genome size
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);
}

sub readbed_list_wps{
	my ($bed_list,$min,$max,$k,$step,$genome,$tgs)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	my %genome_size = %{$tgs};
	open(FILE, $bed_list) || die "BED list $!\n";
	while(chomp(my $bed_file = <FILE>)){
		if($bed_file=~/gz$/){
			open(BED,"gunzip -c $bed_file |") || die "cannot open pipe to $bed_file";
		}else{
			open(BED,$bed_file) || die "BED $!\n";
		}
		while(my $bed_line=<BED>){
			@temp = split /[\ \s\n\t]+/, $bed_line;
			my $len = $temp[2]-$temp[1];
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			#print "$temp[0] GS: $genome_size{$temp[0]}\n";
			if($#temp < 2){
				print STDERR "Not regular BED line?\n$bed_line\n"; 
			}elsif($len>=$min && $len<=$max && exists($genome_size{$temp[0]})){
				$ct++;
				my $lower = int(($temp[1]/$step) - 0.5)*$step;
				my $upper = int(($temp[2]/$step) + 0.5)*$step;
				
				my $left_pos  = $lower + $k;
				my $right_pos = $upper - $k;

				print "$temp[1] $temp[2] $lower $upper $left_pos $right_pos\n";
				for($i=$lower-$k;$i<=$lower+$k;$i+=$step){
					if($i<= int($genome_size{$temp[0]})){
						$twig{$chrid}{$i}--;
					}
				}
				for($i=$upper-$k;$i<=$upper+$k;$i+=$step){
					if($i<= int($genome_size{$temp[0]})){
						$twig{$chrid}{$i}--;
					}
				}
				for($i=$left_pos;$i<=$right_pos;$i+=$step){
					if($i<= int($genome_size{$temp[0]})){
						$twig{$chrid}{$i}++;
					}
				}
			}	
		}
		close(BED);
		print STDERR "Lines used: $ct\n";
	}
	my $GN = 139712364; # Genome size
	if($genome eq 'mm10'){
		$GN = 2800000000;
	}elsif($genome eq 'hg38'){
		$GN = 3300000000;
	}
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}/$ct;
			#my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);

}


sub readbed_list{
	my ($bed_list,$min,$max,$step,$genome,$tgs)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	my %genome_size = %{$tgs};
	open(FILE, $bed_list) || die "BED list $!\n";
	while(chomp(my $bed_file = <FILE>)){
		if($bed_file=~/gz$/){
			open(BED,"gunzip -c $bed_file |") || die "cannot open pipe to $bed_file";
		}else{
			open(BED,$bed_file) || die "BED $!\n";
		}
		while(my $bed_line=<BED>){
			@temp = split /[\ \s\n\t]+/, $bed_line;
			my $len = $temp[2]-$temp[1];
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			#print "$temp[0] GS: $genome_size{$temp[0]}\n";
			if($#temp < 2){
				print STDERR "Not regular BED line?\n$bed_line\n"; 
			}elsif($len>=$min && $len<=$max && exists($genome_size{$temp[0]})){
				my $st=int( ($temp[1] /$step) + 0.5 )*$step;
				my $en=int( ($temp[2] /$step) + 0.5 )*$step;
				for($i=$st;$i<=$en;$i+=$step){
					#print STDERR "$chrid\t$len\t$min\t$max\n";
					if($i<= int($genome_size{$temp[0]})){
						$twig{$chrid}{$i}++;
						$ct++;
					}
				}
			}	
		}
		close(BED);
		print STDERR "Lines used: $ct\n";
	}
	my $GN = 139712364; # Genome size
	if($genome eq 'mm10'){
		$GN = 2800000000;
	}elsif($genome eq 'hg38'){
		$GN = 3300000000;
	}
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);
}

sub readbed_spike{
	my ($bed_file,$min,$max,$step,$genome,$tgs,$spike)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	my %genome_size = %{$tgs};
	open(BED,$bed_file) || die "BED $!\n";
	while(my $bed_line=<BED>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $len = $temp[2]-$temp[1];
		$temp[0]=~s/chr//;
		$chrid=$temp[0];
		#print "$temp[0] GS: $genome_size{$temp[0]}\n";
		if($#temp < 2){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif($len>=$min && $len<=$max && exists($genome_size{$temp[0]})){
			my $st=int( ($temp[1] /$step) + 0.5 )*$step;
			my $en=int( ($temp[2] /$step) + 0.5 )*$step;
			for($i=$st;$i<=$en;$i+=$step){
				#print STDERR "$chrid\t$len\t$min\t$max\n";
				if(exists $genome_size{$temp[0]} && $i<= int($genome_size{$temp[0]})){
					$twig{$chrid}{$i}+= (10000/$spike);
					$ct++;
				}
			}
		}
	}
	close(BED);
	print STDERR "Lines used: $ct\n";
	return (\%twig);
}


sub readbed{
	my ($bed_file,$min,$max,$step,$genome,$tgs)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	my %genome_size = %{$tgs};
	open(BED,$bed_file) || die "BED $!\n";
	while(my $bed_line=<BED>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		my $len = $temp[2]-$temp[1];
		$temp[0]=~s/chr//;
		$chrid=$temp[0];
		#print "$temp[0] GS: $genome_size{$temp[0]}\n";
		if($#temp < 2){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif($len>=$min && $len<=$max && exists($genome_size{$temp[0]})){
			my $st=int( ($temp[1] /$step) + 0.5 )*$step;
			my $en=int( ($temp[2] /$step) + 0.5 )*$step;
			for($i=$st;$i<=$en;$i+=$step){
				#print STDERR "$chrid\t$len\t$min\t$max\t$genome_size{$temp[0]}\n";
				if($i<= int($genome_size{$temp[0]})){
					$twig{$chrid}{$i}++;
					$ct++;
				}
			}
		}
	}
	close(BED);
	print STDERR "Lines used: $ct\n";
	my $GN = 139712364; # Genome size
	if($genome eq 'mm10'){
		$GN = 2800000000;
	}elsif($genome eq 'hg38'){
		$GN = 3300000000;
	}elsif($genome eq 'sacCer3'){
		$GN = 12157105;
	}
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);
}

sub readpairs_fillten{
	my ($pair_file,$min,$max)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0;
	my $i; my $j;
	open(PAIRS,$pair_file) || die "WIG $!\n";
	while(my $bed_line=<PAIRS>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}elsif($temp[5]>=$min && $temp[5]<=$max){
			$temp[0]=~s/chr//;
			$chrid=$temp[0];
			my $lower=int(($temp[2]/10) + 0.5)*10 + 1; 
			my $upper=int(($temp[3]/10) + 0.5)*10 + 1;
			for($i=$lower;$i<=$upper;$i+=10){
				$twig{$chrid}{$i}++;
				$ct++;
			}
		}
	}
	print STDERR "Lines used: $ct\n";
	my $GN = 139712364; # Genome size
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	close(PAIRS);
	return (\%wig);
}

sub pairs2bed{
	# jnk = pairs2bed(pairs_file, out_file (bed) )
	my ($pair_file,$outfile)=@_;
	my @temp; my $bed_line; my $chrid;
	my $i; my $j;
	open(OUT,'>>',$outfile) || die "OUT writeBED $outfile $!\n";
	open(PAIRS,$pair_file) || die "WIG $!\n";
	while(my $bed_line=<PAIRS>){
		@temp = split /[\ \s\n\t]+/, $bed_line;
		if($#temp != 5){
			print STDERR "Not regular BED line?\n$bed_line\n"; 
		}else{
			print OUT "$temp[0]\t$temp[2]\t$temp[3]\n";
		}
	}
	close(PAIRS);
	close(OUT);
	print STDERR "Done with pairs file\n";
	return (1);
}


sub readmultipairs_fillten{
	my ($pair_list,$min,$max)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0; my $read_count=0;
	my $i; my $j;
	open(LIST,$pair_list) || die "PAIR LIST $!\n";
	while(my $pair_file=<LIST>){
		chomp($pair_file);
		print STDERR "opening $pair_file\n";
		open(PAIRS,$pair_file) || die "WIG $!\n";
		while(my $bed_line=<PAIRS>){
			@temp = split /[\ \s\n\t]+/, $bed_line;
			if($#temp != 5){
				print STDERR "Not regular BED line?\n$bed_line\n"; 
			}elsif($temp[5]>=$min && $temp[5]<=$max){
				$read_count++;
				$temp[0]=~s/chr//;
				$chrid=$temp[0];
				my $lower=int(($temp[2]/10) + 0.5)*10 + 1; 
				my $upper=int(($temp[3]/10) + 0.5)*10 + 1;
				for($i=$lower;$i<=$upper;$i+=10){
					$twig{$chrid}{$i}++;
					$ct++;
				}
			}
		}
		close(PAIRS);
		print STDERR "Closing. Read count : $read_count\n";
	}
	close(LIST);
	print STDERR "Positions counted: $ct\n";
	my $GN = 139712364; # Genome size
	#normalize
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			my $normval = $twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);
}
sub readmultipairs_fillone{
	my ($pair_list,$min,$max)=@_;
	my @temp; my %wig; my $bed_line; my $chrid; my %twig; my $ct=0; my $tct=0; my $read_count=0;
	my $i; my $j; my %single_twig;
	open(LIST,$pair_list) || die "PAIR LIST $!\n";
	while(my $pair_file=<LIST>){
		chomp($pair_file);
		print STDERR "opening $pair_file\n";
		open(PAIRS,$pair_file) || die "WIG $!\n";
		while(my $bed_line=<PAIRS>){
			@temp = split /[\ \s\n\t]+/, $bed_line;
			if($#temp != 5){
				print STDERR "Not regular BED line?\n$bed_line\n"; 
			}elsif($temp[5]>=$min && $temp[5]<=$max){
				$read_count++;
				print STDERR "$read_count\n" if($read_count%1000000 == 0);
				$temp[0]=~s/chr//; 
				$chrid=$temp[0];
				my $lower=int(($temp[2]/10) + 1)*10 + 1;
				my $upper=int(($temp[3]/10) - 1)*10 + 1;
				for($i=$temp[2];$i<$lower;$i++){
					$single_twig{$chrid}{$i}++;
				}
				for($i=$upper+1;$i<=$temp[3];$i++){
					$single_twig{$chrid}{$i}++;
				}
				for($i=$lower;$i<=$upper-10;$i+=10){
					$twig{$chrid}{$i}++;
					$tct++;
				}
			}
		}
		close(PAIRS);
		print STDERR "Closing. Read count : $read_count\n";
	}
	close(LIST);
	print STDERR "Positions counted: $ct 10: $tct\n";
	my $GN = 139712364; # Genome size
	#fill up single_twig
	$ct+=(10*$tct);
	foreach $i (keys(%twig)){
		foreach $j ( keys(%{$twig{$i}}) ) {
			for(my $k=$j; $k<=$j+10;$k++){
				$single_twig{$i}{$k}+=$twig{$i}{$j};
			}
		}
	}
	#normalize
	foreach $i (keys(%single_twig)){
		foreach $j ( keys(%{$single_twig{$i}}) ) {
			my $normval = $single_twig{$i}{$j}*$GN/$ct;
			$wig{$i}{$j}=$normval;
		}
	}
	return (\%wig);
}



sub readwig{
#wig format: Coord Value ; chr in the header only
	my @temp; my %wig; my $wig_line; my $chrid;
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//g;
			$chrid=$temp[2];
			print STDERR "$chrid\n";
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $Ywig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	return (\%wig);
}

sub readwig_window{
#wig format: Coord Value ; chr in the header only
	my @temp; my %wig; my $wig_line; my $chrid; my $counter;
	my $window = 10;
	$window = $_[1] if($_[1]);
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//g;
			$chrid=$temp[2];
			print STDERR "$chrid\n";
		}elsif($wig_line=~/^\d/){
			my $pos = (int($temp[0]/$window+0.5))*$window+1;
			$wig{$chrid}{$pos}+=$temp[1]; # $Ywig{Chromosome}{coordinate} = value
		}
	}
	close(WIG);
	return (\%wig);
}


sub readwigArray{
#wig format: Coord Value ; chr in the header only
	my @temp; my %wig; my $wig_line; my $chrid; my @wig_ar;
	open(WIG,$_[0]) || die "WIG $!\n";
	while(my $wig_line=<WIG>){
		@temp = split /[=\ \s\n\t]+/, $wig_line;
		if($wig_line=~/chrom/){
			$temp[2]=~s/chr//;
			$chrid=$temp[2];
		}elsif($wig_line=~/^\d/){
			$wig{$chrid}{$temp[0]}=$temp[1]; # $Ywig{Chromosome}{coordinate} = value
			push(@wig_ar,$temp[1]);
		}
	}
	close(WIG);
	return (\%wig,\@wig_ar);
}

sub readFlank{
#flank format:
# 0 - chr
# 1 - left flank  (genomic, not based on direction of transcription)
# 2 - right flank ( " ")
# 3 - FBgn name
# 4 - filler
# 5 - strand
# 6 - start (TSS)
# 7 - end (TES)
# 8 - flank size (optional)
	my @ftemp; my %gene; my $fline;
	open(FLANKFILE,$_[0]);
	while(chomp($fline=<FLANKFILE>)){
		@ftemp = split /[\ \s\n\t]+/, $fline;
		$gene{$ftemp[0]}{$ftemp[3]}{'start'}  = $ftemp[6];
		$gene{$ftemp[0]}{$ftemp[3]}{'end'}    = $ftemp[7];
		$gene{$ftemp[0]}{$ftemp[3]}{'strand'} = $ftemp[5];
		#convert genomic direction into transcription direction for flanks
		if($ftemp[5] eq "+"){
			$gene{$ftemp[0]}{$ftemp[3]}{'lf'} = $ftemp[1];
			$gene{$ftemp[0]}{$ftemp[3]}{'rf'} = $ftemp[2];
		}
		if($ftemp[5] eq "-"){
			$gene{$ftemp[0]}{$ftemp[3]}{'lf'} = $ftemp[2];
			$gene{$ftemp[0]}{$ftemp[3]}{'rf'} = $ftemp[1];
		}
	}
	close(FLANKFILE);
	return(\%gene);
}

sub normTSS_bp{
# Genehash, wighash, left flank size, right flank size
	my ($href1,$href2,$Lflank,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $lst; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$lst=$thash{$i}{$j}{'start'} + 25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}+$Rflank if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$rst;$m++){
					$count++;
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} +\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$lst=$thash{$i}{$j}{'start'}-25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}-$Rflank if($thash{$i}{$j}{'start'}-$rst > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$rst;$m--){
					my $val = 0;
					$count++;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} -\n";
			}
		}
	}
	return(\%thash);
}

sub normTSS_bp_str_nuc{
# Genehash, nuchash, nwighash, pwighash, left flank size, right flank size
# nuchash --> nucpos{chr}{gn} 
	my ($href1,$href2,$href3,$href4,$Lflank,$Rflank)=@_;
	my %thash = %{$href1};
	my %nuchash = %{$href2};
	my %pMseq  = %{$href3};
	my %nMseq  = %{$href4};
	my $lst; my $rst;
	foreach my $i ( keys(%nuchash) ) {
		foreach my $j ( keys(%{$nuchash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$lst=$thash{$i}{$j}{'start'} + 25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}+$Rflank if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$rst;$m++){
					$count++;
					my $val = 0;
					if(exists $pMseq{$i}{$m}){
						$val=$pMseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				my $str = sprintf("%s\t%s\t+\t%10.4f\t%4d\t",$i,$j,$thash{$i}{$j}{'norm'},$count);
				print STDERR $str;
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$nuchash{$i}{$j}-75;$m++){
					$count++;
					my $val = 0;
					if(exists $pMseq{$i}{$m}){
						$val=$pMseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'st_norm'} = $norm;
				my $str = sprintf("TSS_to_nuc_entry\t%10.4f\t%4d\t",$thash{$i}{$j}{'st_norm'},$count);
				print STDERR $str;
				my $sum=0;
				my $count=0;
				$lst = $nuchash{$i}{$j}-75 if($lst<$nuchash{$i}{$j}-75);
				for(my $m=$lst;$m<=$nuchash{$i}{$j};$m++){
					$count++;
					my $val = 0;
					if(exists $pMseq{$i}{$m}){
						$val=$pMseq{$i}{$m};
					}
					$sum+=$val;
				}
				$norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'nuc_norm'} = $norm;
				my $str = sprintf("TSS_to_nuc_entry\t%10.4f\t%4d\n",$thash{$i}{$j}{'nuc_norm'},$count);
				print STDERR $str;
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$lst=$thash{$i}{$j}{'start'}-25;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}-$Rflank if($thash{$i}{$j}{'start'}-$rst > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$rst;$m--){
					my $val = 0;
					$count++;
					if(exists $nMseq{$i}{$m}){
						$val=$nMseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				my $str = sprintf("%s\t%s\t-\t%10.4f\t%4d\t",$i,$j,$thash{$i}{$j}{'norm'},$count);
				print STDERR $str;
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$nuchash{$i}{$j}+75;$m--){
					my $val = 0;
					$count++;
					if(exists $nMseq{$i}{$m}){
						$val=$nMseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm_st'} = $norm;
				my $str = sprintf("TSS_to_nuc_entry\t%10.4f\t%4d\t",$thash{$i}{$j}{'norm_st'},$count);
				print STDERR $str;
				my $sum=0;
				my $count=0;
				$lst = $nuchash{$i}{$j}+75 if($lst > $nuchash{$i}{$j}+75);
				for(my $m=$lst;$m>=$nuchash{$i}{$j};$m--){
					my $val = 0;
					$count++;
					if(exists $nMseq{$i}{$m}){
						$val=$nMseq{$i}{$m};
					}
					$sum+=$val;
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'nuc_norm'} = $norm;
				my $str = sprintf("TSS_to_nuc_entry\t%10.4f\t%4d\n",$thash{$i}{$j}{'nuc_norm'},$count);
				print STDERR $str;

			}
		}
	}
	return(\%thash);
}


sub normGB_bp{
# Genehash, wighash, Start_Pos, Flank size
	my ($href1,$href2,$Lst,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $st; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$st=$thash{$i}{$j}{'start'} + $Lst;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}+$Rflank if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$st;$m<=$rst;$m++){
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
						$count++;
						$sum+=$val;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} +\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$st=$thash{$i}{$j}{'start'}-$Lst;
				$rst=$thash{$i}{$j}{'rf'};
				$rst=$thash{$i}{$j}{'start'}-$Rflank if($thash{$i}{$j}{'start'}-$rst > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$st;$m>=$rst;$m--){
					my $val = 0;
					if(exists $Mseq{$i}{$m}){
						$val=$Mseq{$i}{$m};
						$count++;
						$sum+=$val;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'} -\n";
			}
		}
	}
	return(\%thash);
}


sub normTSS{
# Genehash, wighash, left flank size, right flank size
	my ($href1,$href2,$Lflank,$Rflank)=@_;
	my %thash = %{$href1};
	my %Mseq  = %{$href2};
	my $lst; my $rst;
	foreach my $i ( keys(%thash) ) {
		foreach my $j ( keys(%{$thash{$i}}) ){
			if($thash{$i}{$j}{'strand'} eq "+"){
				$lst=(int($thash{$i}{$j}{'lf'}/10))*10+1;
				$lst=(int( ($thash{$i}{$j}{'start'}-$Lflank)/10))*10+1 if($thash{$i}{$j}{'start'}-$lst > $Lflank);
				$rst=(int($thash{$i}{$j}{'rf'}/10))*10+1;
				$rst=(int( ($thash{$i}{$j}{'start'}+$Rflank)/10))*10+1 if($rst - $thash{$i}{$j}{'start'} > $Rflank);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m<=$rst;$m+=10){
					if(exists $Mseq{$i}{$m}){
						$sum+=$Mseq{$i}{$m};
						$count++;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'}\n";
			}elsif($thash{$i}{$j}{'strand'} eq "-"){
				$lst=(int($thash{$i}{$j}{'lf'}/10))*10+1;
				$lst=(int( ($thash{$i}{$j}{'start'}+$Lflank)/10))*10+1 if($lst - $thash{$i}{$j}{'start'} > $Lflank);
				$rst=(int($thash{$i}{$j}{'rf'}/10))*10+1;
				$rst=(int( ($thash{$i}{$j}{'start'}-2000)/10))*10+1 if($thash{$i}{$j}{'start'}-$rst > 2000);
				my $sum=0;
				my $count=0;
				for(my $m=$lst;$m>=$rst;$m-=10){
					if(exists $Mseq{$i}{$m}){
						$sum+=$Mseq{$i}{$m};
						$count++;
					}
				}
				my $norm=0;
				$norm=$sum/$count if($count!=0);
				$thash{$i}{$j}{'norm'} = $norm;
				print STDERR "$i $j $thash{$i}{$j}{'norm'}\n";
			}
		}
	}
	return(\%thash);
}

sub median {
  my @ar     = sort { $a <=> $b } @{$_[0]};
  my $mp     = int ( $#ar / 2);
  my $median = $ar[$mp];
	my $num = $#ar+1;
  return ($median,$num);
}



sub mean_se {
  my @ar = @{$_[0]};
  my $sum=0;
	my $x;
  foreach my $ii (@ar){
    $sum+=$ii;
  }
  my $mean=$sum/($#ar+1);
  $sum=0;
  foreach my $ii (@ar){
    $x=($ii-$mean)**2;
    $sum+=$x;
  }
  my $se=(sqrt($sum))/($#ar+1);
	my $num = $#ar+1;
  return ($mean,$se,$num);
}

sub mean_sd {
  my @ar = @{$_[0]};
  my $sum=0;
	my $x;
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


1;
