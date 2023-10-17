#!/usr/bin/perl
#use warnings;

use lib '.';

#note that DPR is not included!
#uncomment the relevant lines in order to calculate and print the DPR
package PrintOut;
my ($outFile, $logFile);

sub FileName
{
	my ($out) = @_;
	$outFile = $out;
}

sub logName
{
	my ($log) = @_;
	$logFile = $log;
}

sub PrintMIME 
{
	print "Content-type:text/html\r\n\r\n";
}

sub PrintInput
{
	my ($out, $seqPrint) = @_;
	if ($seqPrint !~/^>/ && length($seqPrint)>50) { $seqPrint =~ s/(.{50})/$1<br>/g; print $seqPrint, "\n";}
	else {print $seqPrint, "<br>\n"; }
}

sub legend
{
	#print $outFile "chr\tstart\tend\tstrand\tname\tGAGA\tBREu\tTATA\tBREd\tXCPE1\tXCPE2\tMotif1\tdTCT\thTCT\tBBCABW\thInr\tdInr\tMTE\tbridge\tDPE\tDPR\tIIARE\tPB\n"; 
	print $outFile "chr\tstart\tend\tstrand\tname\tGAGA\tBREu\tTATA\tBREd\tXCPE1\tXCPE2\tMotif1\tdTCT\thTCT\tBBCABW\thInr\tdInr\tMTE\tbridge\tDPE\tPB\n"; 
} 

sub legendNoChr
{
	#print $outFile "name\tGAGA\tBREu\tTATA\tBREd\tXCPE1\tXCPE2\tMotif1\tdTCT\thTCT\tBBCABW\thInr\tdInr\tMTE\tbridge\tDPE\tDPR\tIIARE\tPB\n"; 
	print $outFile "name\tGAGA\tBREu\tTATA\tBREd\tXCPE1\tXCPE2\tMotif1\tdTCT\thTCT\tBBCABW\thInr\tdInr\tMTE\tbridge\tDPE\tPB\n"; 
} 

sub input
{
	my ($chr, $start, $end, $strand, $name) = @_;
	print $outFile "$chr\t$start\t$end\t$strand\t$name";
	#print "$chr\t$start\t$end\t$strand\t$name";
}

sub inputNoChr
{
	my ($name) = @_;
	#print $outFile "$chr\t$start\t$end\t$strand\t$name";
	print $outFile  "$name";
}

sub Results
{
	my ($BREu, $TATA, $BREd, $hInr, $dInr, $DPE, $bridge, $MTE, $bridge, $hTCT, $dTCT, $XCPE1, $XCPE2, $PB, $GAGA, $IIARE, $DPR, $BBCABW, $Motif1) = @_; 
	#print $outFile "\t$GAGA\t$BREu\t$TATA\t$BREd\t$XCPE1\t$XCPE2\t$Motif1\t$dTCT\t$hTCT\t$BBCABW\t$hInr\t$dInr\t$MTE\t$bridge\t$DPE\t$DPR\t$IIARE\t$PB\n";
	print $outFile "\t$GAGA\t$BREu\t$TATA\t$BREd\t$XCPE1\t$XCPE2\t$Motif1\t$dTCT\t$hTCT\t$BBCABW\t$hInr\t$dInr\t$MTE\t$bridge\t$DPE\t$PB\n";
}

sub CloseRow
{	
	my ($out) = @_;
	#print "</tr>\n";
	print "\n";
}

sub CloseTable 
{
	my ($out) = @_;
	print "</TABLE>\n";
}

sub NoElements
{
	print "\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0n\to\to"; # 15 times for each element
}

sub NoCombinations
{
	print "<h3> No combinations of core promoter elements were detected in the input sequence. </h3>";
}

sub PrintSequenceTooLong
{
	my ($head, $out) = @_;
	print "\n", "$head is >1000bp\n";
	print "Please use sequence length of less than 1000bp.\n";
	
	print $logFile "\nseq >1000bp: $head\n";
}

sub PrintOutOfScore
{
	my ($out, $name, $cutoff, $lower) = @_;
	
	if ($cutoff < $lower or $cutoff > 20)  #GC_04_2023 - update upper limit verificattio of cutoff 2 supporl log values
	{
		print "\n", "The ", $name, " score you have entered is outside the allowed range!\n";
		print "Please use a value between ", $lower, " and 1 \n";
		print $logFile "Please use a cutoff value between 0 and 1 for $name\n";
		return 1;
		die;
	}
	else
	{
		return 0;
	}
}

sub PrintOutOfGCrange # GC_04_2023
{
	my ($gc) = @_;
	if ($gc< 1 or $gc > 100)  #check whether this is reasonale range
	{
		print "\n", "The value of GC content ", $gc, " you have entered is outside the allowed range!\n";
		print "Please use a value between 10 and 90 \n\n ";
		return 1;
		#die;

	}
	else { return 0}
}
sub CharWarning
{
	my ($out, $char, $CurSeq, $head) = @_;
	print "\nThe use of characters other than ACGT was detected!\n";
	#print "Please verify your input sequence contains ACGT characters only.\n";
	print "The problematic character is ", $char, "\n";
	print "within the following sequence \n$head\n$CurSeq\n";
	#print "within the following sequence (bare)\n", $CurSeq;
}

sub ProblemChar
{
	my ($out, $head, $char) = @_;

	print "non ACGT character: ", $char, "\n";
	print "Problematic seq: $head\n";
	
	print $logFile "\nThe use of characters other than ACGT was detected!\n";
	print $logFile "Problematic character: ", $char, "\n";
	print $logFile "Problematic seq: $head\n";
}

sub Cite
{
	print $logFile "\nPlease cite the following paper if you found this program useful- \n";
	print $logFile "ElemeNT: A Computational Tool for Detecting Core Promoter Elements.\n";
	print $logFile "Sloutskin A, Danino YM, Orenstein Y, Zehavi Y, Doniger T, Shamir R, and Juven-Gershon T.\n";
	print $logFile "Transcription (2015)\n";
	
	print "ElemeNT done\n"; #command line print for QC
}

sub logStart
{
	my ($curTime,$mday,$mon,$year,$day) = @_;
	print $logFile "ElemeNT command line tool V2023\n";  #GC_04_2023
	print $logFile "local date and time: $day, $mday $mon $year, $curTime\n";
	#print $logFile "\nNote: the number and names of the processed sequences will be displayed at the end if no fatal error is encountered\n"; 
}

sub logParam
{
	my ($inputName, $outputName, $configStr, $chrUse, $smoothCon, $startCon, $revCom, $gc, @Order) = @_;
	print $logFile "\ninput: $inputName\n";
	print $logFile "output: $outputName\n";
	print $logFile "run on anti-sense seq: $revCom\n";
	print $logFile "extract coordinates: $chrUse\n";
	print $logFile "TSS at position: $startCon\n";
	print $logFile "smoothing window (bp): $smoothCon\n";
	print $logFile "GC content (%): $gc\n";
	
	my $dep = 0;
	
	print $logFile "\nprocessed elements and cutoff scores\n";
	foreach (@Order)	#printing just the relevant elements and cutoffs
	{
		$name = $_;
		$searchName = $name."Cutoff";
		#$configStr =~ /$searchName: (\d*\.?\d*);/;
		$configStr =~ /$searchName: (-?[0-9]\d*(\.\d+));/;  #GC_04_2023
		print $logFile "$name: $1\n";

		if ($name eq "hInr") #printing hInr dependencies
		{
			$dep++;
			my @Dep = ("MTE", "bridge1", "bridge2", "DPE");
			foreach (@Dep)
			{
				$name = $_;
				$searchName = $name."Cutoff";
				if ($configStr =~ /$searchName: (\d*\.?\d*);/) {print $logFile "$name: $1\n";}

			}
		}
		if ($dep == 0 && $name eq "dInr") #printing dInr dependencies in case no hInr is specified 
		{
			my @Dep = ("MTE", "bridge1", "bridge2", "DPE");
			foreach (@Dep)
			{
				$name = $_;
				$searchName = $name."Cutoff";
				if ($configStr =~ /$searchName: (\d*\.?\d*);/) {print $logFile "$name: $1\n";}

			}
		}
	}
	print $logFile "\n";
}

sub logNum
{
	my ($numSeq) = @_;
	print $logFile "\n$numSeq sequences processed\n";
}

sub logSeq 
{
	my ($name) = @_;
	print $logFile "processed successfully: $name\n";
}

sub logNoElem
{
	my ($elem, $cwd) = @_;
	print $logFile "\ncan't find $elem.txt in folder\n$cwd/Elements";
}

1;
