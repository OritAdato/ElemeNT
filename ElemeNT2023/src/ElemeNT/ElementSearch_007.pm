#!/usr/bin/perl
#use warnings;

use lib '.';
use Cwd;
my $cwd = getcwd();
use PrintOut_007;

package ElemeNTSearch;

$smooth, $start;

my $startOfString = 0;
my $endOfString = 1;

$hInr = 1;
$dInr = 2;
$MTE = 3;
$DPE = 4;
$bridge1 = 5;
$bridge2 = 6;
$hTCT = 7;
$dTCT = 8;
$BREu = 9;
$BREd = 10;
$TATA = 11;
$XCPE1 = 12;
$XCPE2 = 13;
$PB = 14;
$GAGA = 15; #21.4.20
#$IIARE = 16;
#$DPR = 17;
$BBCABW = 18;
$Motif1 = 19;

$lastBridge1Score = 1;

@hInrRes = ();
@dInrRes = ();
@TATARes = ();
@BREuRes = ();
@BREdRes = ();
@hTCTRes = ();
@dTCTRes = ();
@XCPE1Res = ();
@XCPE2Res = ();
@MTERes = ();
@PBRes = (); #5.4.20
@GAGARes = ();
#@IIARERes = ();
#@DPRRes = ();
@BBCABWRes = ();
@Motif1Res = ();

sub ElemeNT
{
	#maybe take out $freqStr
	my ($head, $out, $fh, $elementType, $freqStr, $configStr, $gc, $InitiatorPos, $ii, $InrType) = @_;
	my ($elementName, $cutoff, $lower, $upper, $InitiatorPosOffset) = Init ($elementType, $configStr);

	if (PrintOut::PrintOutOfScore ($out, $elementName, $cutoff, $lower) == 0)
	{
#		if ($InitiatorPos + $InitiatorPosOffset > length ($fh)) {next;}
		if ($InitiatorPos + $InitiatorPosOffset <= length ($fh)) {
			my $patternStr = GetPatternStr ($elementType, $fh, $InitiatorPos + $InitiatorPosOffset, $upper);
			my @ResOut;
			my @ScoreRes = GetScore ($elementType, $patternStr, $freqStr, $cutoff, $out, $head, $gc);
			if ($elementType == $bridge1) { $lastBridge1Score = $ScoreRes[0];}
			else {	@ResOut = SortRes ($cutoff, @ScoreRes); }

			if ($elementType == $hInr) { @hInrRes = RememberArray(@ResOut); }
			elsif ($elementType == $dInr) { @dInrRes = RememberArray(@ResOut);}
			elsif ($elementType == $TATA) { @TATARes = RememberArray(@ResOut);}
			elsif ($elementType == $BREu) { @BREuRes = RememberArray(@ResOut);}
			elsif ($elementType == $BREd) { @BREdRes = RememberArray(@ResOut);}
			elsif ($elementType == $hTCT) { @hTCTRes = RememberArray(@ResOut);}
			elsif ($elementType == $dTCT) { @dTCTRes = RememberArray(@ResOut);}
			elsif ($elementType == $XCPE1) { @XCPE1Res = RememberArray(@ResOut);}
			elsif ($elementType == $XCPE2) { @XCPE2Res = RememberArray(@ResOut);}
			elsif ($elementType == $PB) { @PBRes = RememberArray(@ResOut);} #5.4.20
			elsif ($elementType == $GAGA) { @GAGARes = RememberArray(@ResOut);}
		#elsif ($elementType == $IIARE) { @IIARERes = RememberArray(@ResOut);}
		#elsif ($elementType == $DPR) { @DPRRes = RememberArray(@ResOut);}
			elsif ($elementType == $BBCABW) { @BBCABWRes = RememberArray(@ResOut);}
			elsif ($elementType == $Motif1) { @Motif1Res = RememberArray(@ResOut);}

			if ($elementType == $MTE || $elementType == $DPE || $elementType == $bridge2)
			{
				RemDep($elementType, $InrType, $ResOut[0][0], $InitiatorPosOffset, $ii)
			}

		}
	}
}
#input: element type and configStr. process and output: calculates the cutoff according to the element type and returns the element name, cutoff and 3 numeric values (the 1st represents the distance in relation to the TSS)
sub Init
{
	my ($elementType, $configStr) = @_;
	my ($cutoff) = GetConfigCutoff($configStr, $elementType);
	if ($elementType == $hInr) {return "Mammalian Initiator", $cutoff, -10.2, 7, 0;}
	if ($elementType == $dInr) {return "Drosophila Initiator", $cutoff, -14.7, 6, 0;}
	if ($elementType == $MTE)  {return "MTE", $cutoff, -34.6, 12, 19;}
	if ($elementType == $DPE)  {return "DPE", $cutoff, -14.9, 6, 29;}
	if ($elementType == $bridge1)  { return "bridge 1", $cutoff, -3.2, 5, 19;}
	if ($elementType == $bridge2)  { return "bridge 2", $cutoff, -12.7, 4, 19+12;}
	if ($elementType == $hTCT) {return "Human TCT", $cutoff, -23.7, 7, 0;}
	if ($elementType == $dTCT) {return "Drosophila TCT", $cutoff, -30.9, 8, 0;}
	if ($elementType == $BREu) {return "BRE upstream", $cutoff, -27.2, 7, 0;}
	if ($elementType == $BREd) {return "BRE downstream", $cutoff, -9.4, 7, 0;}
	if ($elementType == $TATA) {return "TATA box", $cutoff, 2.7, 8, 0;}
	if ($elementType == $XCPE1) {return "XCPE 1", $cutoff, -48, 10, 0;}
	if ($elementType == $XCPE2) {return "XCPE 2", $cutoff, -42.7, 11, 0;}
	if ($elementType == $PB) {return "Pause button", $cutoff, -40.2, 7, 0;} #5.4.20
	if ($elementType == $GAGA) {return "GAGA factor", $cutoff, -41.3, 10, 0;}
	#if ($elementType == $IIARE) {return "TFIIA response element", $cutoff, 0, 10, 0;}
	#if ($elementType == $DPR) {return "DPR element", $cutoff, 0, 19, 0;}
	if ($elementType == $BBCABW) {return "BBCABW Initiator", $cutoff, -11.5, 6, 0;}
	if ($elementType == $Motif1) {return "Motif1 element", $cutoff, -100.4, 14, 0;}
}

#input: configStr and element type (a number). process and output: defines the cutoff type according to the element type and returns the cutoff value that is in the configStr
sub GetConfigCutoff
{
	my ($configStr, $elementType) = @_;
	my $searchName = "";

	if ($elementType == $hInr) { $searchName = "hInrCutoff"; }
	if ($elementType == $dInr) { $searchName = "dInrCutoff"; }
	if ($elementType == $MTE) { $searchName = "MTECutoff"; }
	if ($elementType == $DPE) { $searchName = "DPECutoff"; }
	if ($elementType == $bridge1) { $searchName = "bridge1Cutoff"; }
	if ($elementType == $bridge2) { $searchName = "bridge2Cutoff"; }
	if ($elementType == $hTCT) { $searchName = "hTCTCutoff"; }
	if ($elementType == $dTCT) { $searchName = "dTCTCutoff"; }
	if ($elementType == $BREu) { $searchName = "BREuCutoff"; }
	if ($elementType == $BREd) { $searchName = "BREdCutoff"; }
	if ($elementType == $TATA) { $searchName = "TATACutoff"; }
	if ($elementType == $XCPE1) { $searchName = "XCPE1Cutoff"; }
	if ($elementType == $XCPE2) { $searchName = "XCPE2Cutoff"; }
	if ($elementType == $PB) { $searchName = "PBCutoff"; } #5.4.20
	if ($elementType == $GAGA) { $searchName = "GAGACutoff"; }
	#if ($elementType == $IIARE) { $searchName = "IIARECutoff"; }
	#if ($elementType == $DPR) { $searchName = "DPRCutoff"; }
	if ($elementType == $BBCABW) { $searchName = "BBCABWCutoff"; }
	if ($elementType == $Motif1) { $searchName = "Motif1Cutoff"; }

	if ($searchName ne "" && $configStr =~ /$searchName: (-?[0-9]\d*(\.\d+));/) # GC_04_2023
	{
		return ($1);
	}
	return 0;
}

# input: None. process and output: defines and opens the Position Weight Matrix (PWM) according to the element type
sub GetScore
{
	my $skip = 0;
	my ($elementType, $fh, $freqStr, $cutoff, $out, $head, $gc) = @_;
	my $score = 0;
	my $matrixName;

	#my @FreqArray = split ("\t", $freqStr);
	if ($elementType == $dInr) { $matrixName = "dInr";}
	if ($elementType == $hInr) { $matrixName = "hInr";}
	if ($elementType == $TATA) { $matrixName = "TATA";}
	if ($elementType == $hTCT) { $matrixName = "hTCT";}
	if ($elementType == $dTCT) { $matrixName = "dTCT";}
	if ($elementType == $BREu) { $matrixName = "BREu";}
	if ($elementType == $BREd) { $matrixName = "BREd";}
	if ($elementType == $XCPE1) { $matrixName = "XCPE1";}
	if ($elementType == $XCPE2) { $matrixName = "XCPE2";}
	if ($elementType == $MTE) { $matrixName = "MTE";}
	if ($elementType == $DPE) { $matrixName = "DPE";}
	if ($elementType == $bridge1) { $matrixName = "bridge1";}
	if ($elementType == $bridge2) { $matrixName = "bridge2";}
	if ($elementType == $PB) { $matrixName = "PB";} #5.4.20
	if ($elementType == $GAGA) { $matrixName = "GAGA";}
	#if ($elementType == $IIARE) { $matrixName = "IIARE";}
	#if ($elementType == $DPR) { $matrixName = "DPR";}
	if ($elementType == $BBCABW) { $matrixName = "BBCABW";}
	if ($elementType == $Motif1) { $matrixName = "Motif1";}

	if (!open my $matrix, "<", "$cwd/Elements/$matrixName.txt") {
		PrintOut::logNoElem ($matrixName, $cwd);
	}

	open (my $matrix, "<", "$cwd/Elements/$matrixName.txt")
		or die "Can't open $matrixName.txt";
	
	my @lines = <$matrix>;
	my @CellsArray;
	my @TheMatrix;

	#creating the array from the input matrix
	for (my $lnum = 0; $lnum < 4; $lnum++)
	{
		chomp $lines[$lnum];
		@CellsArray = split ("\t", $lines[$lnum]);

		for (my $colnum = 0; $colnum <= $#CellsArray; $colnum++)
		{
			$TheMatrix[$lnum][$colnum] = $CellsArray[$colnum];
		}
	}
	
		# creating PWM
		# for each value, divide by the background frequency and take log2 of the product
	# calculate background frequency
	if (PrintOut::PrintOutOfGCrange ($gc) != 0) {die}
	my $gc_freq = ($gc / 100) / 2;     # GC_04_2023 - freq of g & c
	my $at_freq = (1 - $gc / 100) / 2; # GC_04_2023 -  freq of a & t

	my $bckgrnd_freq = 0.25;
	for (my $lnum = 0; $lnum < 4; $lnum++) 
	{
		if (($lnum == 0) or ($lnum == 3)) {      #GC_04_2023
			$bckgrnd_freq = $at_freq;            #GC_04_2023
		}
		else {
			$bckgrnd_freq = $gc_freq;            #GC_04_2023
		}
		for (my $colnum = 0; $colnum <= $#CellsArray; $colnum++)
		{
			$TheMatrix[$lnum][$colnum] = log2($TheMatrix[$lnum][$colnum]/$bckgrnd_freq); #GC_04_2023 gc background
		}
	}
	
	# #printing the created PWM
	# for my $lnum (@TheMatrix) {  print join(",", @{$lnum}), "\n"; }
	# print "\nenough\n";
	
	# #calculating and the maximum score and normalizing all values to it.
	# my $max_col;

	# for (my $colnum = 0; $colnum <= $#CellsArray; $colnum++)
	# {
		# $max_col = 0;
		# for (my $lnum = 0; $lnum < 4; $lnum++)
		# {
			# if ($TheMatrix[$lnum][$colnum] > $max_col) { $max_col = $TheMatrix[$lnum][$colnum] }
		# }
		# for (my $lnum = 0; $lnum < 4; $lnum++)
		# {
			# $TheMatrix[$lnum][$colnum] = $TheMatrix[$lnum][$colnum]/$max_col;
		# }
	# }
	# calculating the product of probability for each position
	my $basenm;
	my $pos_sh;
	my $Prob_i;
	my $strpos;
	my @ScoreRes;

	for ($strpos = 0; $strpos <= (length($fh) - ($#CellsArray + 1)); $strpos++)
	{
		#$Prob_i = 1.0; GC_04_2023
		$Prob_i = 0;  # GC_04_2023
		for (my $pos_sh = 0; $pos_sh <= $#CellsArray; $pos_sh++)
		{
			if (substr($fh, $strpos + $pos_sh, 1) =~ /A/i) { $basenm = 0; }
			elsif (substr($fh, $strpos + $pos_sh, 1) =~ /C/i) { $basenm = 1; }
			elsif (substr($fh, $strpos + $pos_sh, 1) =~ /G/i) { $basenm = 2; }
			elsif (substr($fh, $strpos + $pos_sh, 1)  =~ /T/i) { $basenm = 3; }
			else #not relevant, checking for ACGT characters before run, 4/5/2020
			{
				PrintOut::CharWarning ($out, substr($fh, $strpos + $pos_sh, 1), $fh, $head);
				die;
			}
		#$Prob_i = $Prob_i * $TheMatrix[$basenm][$pos_sh]; #for probability matrix
		$Prob_i = $Prob_i + $TheMatrix[$basenm][$pos_sh]; #add, for log-based ,matrix
		}

		$ScoreRes[$strpos] = $Prob_i;

		if ($elementType == $bridge2)
		{
			$ScoreRes[$strpos]+= $lastBridge1Score;
		}
	}
	close $matrix;   #GC_04_2023
	return @ScoreRes;

	
	sub log2 
	{
		my $n = shift;
		if ($n == 0) { return -30; } #deal with zero values, assign verly low number
		return log($n)/log(2);
	}
	
}

sub SortRes
{
	my ($cutoff, @ResIn) = @_;
	my @ResOut;
	my $posResIn = 0;
	my $posResOut = 0;
	foreach (@ResIn)
	{
		if ($_ >= $cutoff)
		{
			$ResOut[$posResOut][0] = $_;
			$ResOut[$posResOut][1] = $posResIn;
			$posResOut++;
		}
		$posResIn++;
	}
	if ($#ResOut > 0)
	{
		use sort 'stable';
		@ResOut = reverse sort { $a->[0] <=> $b->[0] } (@ResOut);
	}
	return @ResOut;
}

sub GetDependencies
{
	my ($configStr, $fh, $head, $gc) = @_;

	for (my $ii = 0; $ii <= $#hInrRes; $ii++)
		{
			$hInrRes[$ii][2] = -1; $hInrRes[$ii][3] = -1; #MTE columns
			$hInrRes[$ii][4] = -1; $hInrRes[$ii][5] = -1; #DPE columns
			$hInrRes[$ii][6] = -1; $hInrRes[$ii][7] = -1; #bridge columns, based on br2 score and br1 position

			my $searchName = "";
			if ($configStr =~ /MTECutoff/)
			{
				ElemeNT ($head, $out, $fh, $MTE, $freqStr, $configStr, $gc, $hInrRes[$ii][1], $ii, $hInr);
			}
			if ($configStr =~ /DPECutoff/)
			{
				ElemeNT ($head, $out, $fh, $DPE, $freqStr, $configStr, $gc, $hInrRes[$ii][1], $ii, $hInr);
			}
			if ($configStr =~ /bridge1Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge1, $freqStr, $configStr, $gc, $hInrRes[$ii][1], $ii, $hInr)
			}
			if ($configStr =~ /bridge2Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge2, $freqStr, $configStr, $gc, $hInrRes[$ii][1], $ii, $hInr)
			}
		}

	for (my $ii = 0; $ii <= $#dInrRes; $ii++)
		{
			$dInrRes[$ii][2] = 0; $dInrRes[$ii][3] = 0; #MTE columns
			$dInrRes[$ii][4] = 0; $dInrRes[$ii][5] = 0; #DPE columns
			$dInrRes[$ii][6] = 0; $dInrRes[$ii][7] = 0; #bridge columns, based on br2 score and br1 position

			my $searchName = "";
			if ($configStr =~ /MTECutoff/)
			{
				ElemeNT ($head, $out, $fh, $MTE, $freqStr, $configStr, $gc, $dInrRes[$ii][1], $ii, $dInr);
			}
			if ($configStr =~ /DPECutoff/)
			{
				ElemeNT ($head, $out, $fh, $DPE, $freqStr, $configStr, $gc, $dInrRes[$ii][1], $ii, $dInr);
			}
			if ($configStr =~ /bridge1Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge1, $freqStr, $configStr, $gc, $dInrRes[$ii][1], $ii, $dInr)
			}
			if ($configStr =~ /bridge2Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge2, $freqStr, $configStr, $gc, $dInrRes[$ii][1], $ii, $dInr)
			}
		}
		
		for (my $ii = 0; $ii <= $#BBCABWRes; $ii++) #looking for dependencies from a BBACBw Inr
		{
			$BBCABWRes[$ii][2] = 0; $BBCABWRes[$ii][3] = 0; #MTE columns
			$BBCABWRes[$ii][4] = 0; $BBCABWRes[$ii][5] = 0; #DPE columns
			$BBCABWRes[$ii][6] = 0; $BBCABWRes[$ii][7] = 0; #bridge columns, based on br2 score and br1 position
			
			#updating the position by 1, to account for the 1bp difference between d/hInr and bbCABW (2/3nt before the A+1, respectively)
			my $pos = $BBCABWRes[$ii][1] + 1;  

			my $searchName = "";
			if ($configStr =~ /MTECutoff/)
			{
				ElemeNT ($head, $out, $fh, $MTE, $freqStr, $configStr, $gc, $pos, $ii, $BBCABW);
			}
			if ($configStr =~ /DPECutoff/)
			{
				ElemeNT ($head, $out, $fh, $DPE, $freqStr, $configStr, $gc, $pos, $ii, $BBCABW);
			}
			if ($configStr =~ /bridge1Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge1, $freqStr, $configStr, $gc, $pos, $ii, $BBCABW)
			}
			if ($configStr =~ /bridge2Cutoff/)
			{
				ElemeNT ($head, $out, $fh, $bridge2, $freqStr, $configStr, $gc, $pos, $ii, $BBCABW)
			}
		}
}

sub GetPatternStr
{
	my ($elementType, $fh, $PatternPos, $len) = @_;

	if ($elementType == $MTE || $elementType == $DPE || $elementType == $bridge1 || $elementType == $bridge2)
	{
		return substr ($fh, $PatternPos, $len);
	}
	return $fh;
}

sub InsideOfString
{
	my ($posString, $stringLength, $start) = @_;
	my $isInside = 1;
	if ($start == 1 && $posString < 0)
	{
		$isInside = 0;
	}
	if ($start == 0 && $posString >= $stringLength)
	{
		$isInside = 0;
	}
	return $isInside;
}

sub RememberArray
{
	my (@ResOut) = @_;
	my ($row, $col);

	my @Remember;
	for($row = 0; $row <= $#ResOut; $row++)
	{
		for($col = 0; $col < 2; $col++) { $Remember[$row][$col] = $ResOut[$row][$col]; }
	}
	return @Remember;
}

sub RemDep
{
	my ($elementType, $InrType, $RemScore, $InitiatorPosOffset, $ii) = @_;
	my $InrName;

	if ($InrType == $hInr) {$InrName = "hInrRes";}
	if ($InrType == $dInr) {$InrName = "dInrRes";}
	if ($InrType == $BBCABW) {$InrName = "BBCABWRes"; }

	if ($elementType == $DPE )
	{
		if (defined($RemScore))
		{
			$$InrName[$ii][4] = $RemScore;
			$$InrName[$ii][5] = $$InrName[$ii][1] + $InitiatorPosOffset;
		}
	}
	if ($elementType == $MTE )
	{
		if (defined($RemScore))
		{
			$$InrName[$ii][2] = $RemScore;
			$$InrName[$ii][3] = $$InrName[$ii][1] + $InitiatorPosOffset;
		}
	}
	if ($elementType == $bridge2 )
	{
		if (defined($RemScore))
		{
			$$InrName[$ii][6] = $RemScore;
			$$InrName[$ii][7] = $$InrName[$ii][1] + 19;
		}
	}
}

 sub PrintTable
{
	my ($fh, @Order) = @_;
	my ($offset, $found);
	my ($BREuConc, $TATAConc, $BREdConc, $hInrConc, $dInrConc, $MTEConc, $DPEConc, $bridgeConc, $hTCTConc, $dTCTConc, $XCPE1conc, $XCPE2conc, $PBconc, $GAGAconc, $IIAREconc, $DPRconc, $BBCABWconc, $Motif1conc); #initializing the conclusion values (0 or 1)

	foreach (@Order)
	{
		$found = "";
		if ($_ eq 'hInr') {$offset = -2; $arrName = "hInrRes"; }
		if ($_ eq 'dInr') {$offset = -2; $arrName = "dInrRes"; }
		elsif ($_ eq 'TATA') {$offset = -30; $arrName = "TATARes"; }
		elsif ($_ eq 'hTCT') {$offset = -1; $arrName = "hTCTRes"; }
		elsif ($_ eq 'dTCT') {$offset = -2; $arrName = "dTCTRes"; }
		elsif ($_ eq 'BREu') {$offset = -37; $arrName = "BREuRes"; }
		elsif ($_ eq 'BREd') {$offset = -24; $arrName = "BREdRes"; }
		elsif ($_ eq 'XCPE1') {$offset = -8; $arrName = "XCPE1Res"; }
		elsif ($_ eq 'XCPE2') {$offset = -9; $arrName = "XCPE2Res"; }
		elsif ($_ eq 'PB') {$offset = 25; $arrName = "PBRes";} #5.4.20
		elsif ($_ eq 'GAGA') {$offset = -80; $arrName = "GAGARes";}
		#elsif ($_ eq 'IIARE') {$offset = 24; $arrName = "IIARERes";}
		#elsif ($_ eq 'DPR') {$offset = 17; $arrName = "DPRRes";}
		elsif ($_ eq 'BBCABW') {$offset = -3; $arrName = "BBCABWRes";}
		elsif ($_ eq 'Motif1') {$offset = -7; $arrName = "Motif1Res";}

		#print("arr: ",$arrName," ;");  #GC_04_2023
		for (my $ii = 0; $ii <= $#$arrName; $ii++)
		{
			my $StartPos = $$arrName[$ii][1]; my $score = sprintf "%.2f", $$arrName[$ii][0]; #GC_04_2023 from %.4f-> %.2f
			if ((($start + $offset - $smooth) <= ($StartPos)) &&  (($StartPos) <= ($start + $offset + $smooth)))
			#the position should be within the smoothed range around the pos relative to the start point

			{
				$StartPos -= $start; # to get the relative position
				if ($StartPos >= 0) {$StartPos++; } #to shift the 0-to-positive side to start from +1
				$found .= "$StartPos,$score;";
				
				if ($arrName eq "hInrRes" || $arrName eq "dInrRes" || $arrName eq "BBCABWRes") #printing initiator dependencies
				{
					#assigning Inr type
					my $InrType = "";
					if ($arrName eq "hInrRes")  {$InrType = "hInr";}
					elsif ($arrName eq "dInrRes")  {$InrType = "dInr";}
					elsif ($arrName eq "BBCABWRes")  {$InrType = "BBCABW";}
					
					if ($$arrName[$ii][2] > 0) #MTE
					{
						$score = sprintf "%.2f",$$arrName[$ii][2]; #GC_04_2023 no need to update the start position again
						$MTEConc .= "$StartPos,$score($InrType);";
					}
					if ($$arrName[$ii][4] > 0) #DPE
					{
						$score = sprintf "%.2f",$$arrName[$ii][4]; #GC_04_2023

						$DPEConc .= "$StartPos,$score($InrType);";
					}
					if ($$arrName[$ii][6] > 0) #bridge
					{
						$score = sprintf "%.2f",$$arrName[$ii][6];  #GC_04_2023

						$bridgeConc .= "$StartPos,$score($InrType);";
					}
				}
			}
		}

		my $missing = "no";
		if ($_ eq 'hInr') {$hInrConc = $found; if ($hInrConc eq "") {$hInrConc = $missing;}}
		if ($_ eq 'dInr') {$dInrConc = $found; if ($dInrConc eq "") {$dInrConc = $missing;}}
		if ($_ eq 'TATA') {$TATAConc = $found; if ($TATAConc eq "") {$TATAConc = $missing;}}
		if ($_ eq 'hTCT') {$hTCTConc = $found; if ($hTCTConc eq "") {$hTCTConc = $missing;}}
		if ($_ eq 'dTCT') {$dTCTConc = $found; if ($dTCTConc eq "") {$dTCTConc = $missing;}}
		if ($_ eq 'BREu') {$BREuConc = $found; if ($BREuConc eq "") {$BREuConc = $missing;}}
		if ($_ eq 'BREd') {$BREdConc = $found; if ($BREdConc eq "") {$BREdConc = $missing;}}
		if ($_ eq 'XCPE1') {$XCPE1Conc = $found; if ($XCPE1Conc eq "") {$XCPE1Conc = $missing;}}
		if ($_ eq 'XCPE2') {$XCPE2Conc = $found; if ($XCPE2Conc eq "") {$XCPE2Conc = $missing;}}
		if ($_ eq 'PB') {$PBConc = $found; if ($PBConc eq "") {$PBConc = $missing;}} #5.4.20
		if ($_ eq 'GAGA') {$GAGAConc = $found; if ($GAGAConc eq "") {$GAGAConc = $missing;}}
		#if ($_ eq 'IIARE') {$IIAREConc = $found; if ($IIAREConc eq "") {$IIAREConc = $missing;}}
		#if ($_ eq 'DPR') {$DPRConc = $found; if ($DPRConc eq "") {$DPRConc = $missing;}}
		if ($_ eq 'BBCABW') {$BBCABWConc = $found; if ($BBCABWConc eq "") {$BBCABWConc = $missing;}}
		if ($_ eq 'Motif1') {$Motif1Conc = $found; if ($Motif1Conc eq "") {$Motif1Conc = $missing;}}

		if ($DPEConc =~ /(^$)|($missing)/) {$DPEConc = $missing;}
		if ($MTEConc =~ /(^$)|($missing)/) {$MTEConc = $missing;}
		if ($bridgeConc =~ /(^$)|($missing)/) {$bridgeConc = $missing;}
	}
	PrintOut::Results($BREuConc, $TATAConc, $BREdConc, $hInrConc, $dInrConc, $DPEConc, $bridgeConc, $MTEConc, $bridgeConc, $hTCTConc, $dTCTConc, $XCPE1Conc, $XCPE2Conc, $PBConc, $GAGAConc, $IIAREConc, $DPRConc, $BBCABWConc, $Motif1Conc); #5.4.20- PBConc Addition
}

1;