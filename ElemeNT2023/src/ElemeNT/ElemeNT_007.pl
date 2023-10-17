#!/usr/bin/perl
#use strict;
#use warnings;

use lib '.';
use Cwd();
use ElementSearch_007;
use PrintOut_007;
#getting and printing the local time
my ($monHum, $dayHum, $curTime, $curTimeLog);
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime;
$year += 1900; #updating the year from the stored variable
$monHum = $months[$mon]; #creating human-readable output
$dayHum = $days[$wday]; #creating human-readable output
#print "$mday$months[$mon]$year\_$hour$min$sec\n";

$curTime = sprintf("%02d%02d%02d", $hour, $min, $sec);
$curTimeLog = sprintf("%02d:%02d:%02d", $hour, $min, $sec);

#open (my $log, ">", "ElemeNT.log") #creating a global log file
open (my $log, ">", "ElemeNT\_$mday$monHum$year\_$curTime.log") #creating log file for each run
	   or die "Can't open log_ElemeNT.log";

PrintOut::logName ($log);
PrintOut::logStart ($curTimeLog,$mday,$monHum,$year,$dayHum);
my $cwd = getcwd();     #GC_04_2023
my $num = 0; #initializing sequence counter
open (my $fh1, "<", "Config.txt")
	   or die "Can't open Config.txt";

my ($inputName, $outputName, $organism, $order, $configStr, $chrUse, $smoothCon, $startCon, $revCom, $gc) = GetConfig($fh1);
my @OrderArray = split (", ", $order);

PrintOut::logParam ($inputName, $outputName, $configStr, $chrUse, $smoothCon, $startCon, $revCom, $gc, @OrderArray);

$ElemeNTSearch::smooth = $smoothCon;
$ElemeNTSearch::start = $startCon;

open (my $fh, "<", $inputName)
	   or die "Can't open $inputName";

open (my $out1, ">", $outputName)
	   or die "Can't open $outputName";

PrintOut::FileName ($out1);

sub GetConfig
{
	my $configStr = "";
	my ($inputName, $outputName, $order, $chrUse, $smoothCon, $startCon, $revCom, $gc);
	while (<$fh1>)
	{
		if ($_ =~ /InputFileName:\s(.+);/)  { $inputName = $1; }
		elsif  ($_ =~ /OutputFileName:\s(.+);/) { $outputName = $1; }
		elsif  ($_ =~ /chrData:\s(.+);/) { $chrUse = $1; } #extract UCSC genomic data from name
		elsif  ($_ =~ /ReverseSeq:\s(.+);/) { $revCom = $1; } #whether the sequence should be run also as revcom
		#elsif  ($_ =~ /UseTSS:\s(.+);/) { $useTSS = $1; } #for future use- define if want to use the TSS or not
		elsif  ($_ =~ /SmoothConstant:\s(.+);/) { $smoothCon = $1; } #smooth constant
		elsif  ($_ =~ /StartConstant:\s(.+);/) { $startCon = $1; } #TSS to use
		elsif  ($_ =~ /Organism:\s(.+);/)   { $organism = $1; } #not relevant apr 2020
		elsif ($_ =~ /Order:\s(.+);/)  { $order = $1; }
		elsif ($_ =~ /GC:\s(.+);/)  { $gc = $1; } #GC_04_2023
		else   { chomp; $configStr .= $_ ; }
	}
	return ($inputName, $outputName, $organism, $order, $configStr, $chrUse, $smoothCon, $startCon, $revCom, $gc);
}
	   
my $freqStr = 0;
my $headStr = "";


my $seq = "";
my ($rev, $head);

if ($chrUse eq "yes") { PrintOut::legend; }
elsif ($chrUse eq "no") { PrintOut::legendNoChr; }


while (<$fh>)
{
	chomp $_;
	#using only multiple fasta format- apr2020
	if (/(^>)/) 
	{	
		if (length $seq != 0) 
		{ 
			ElemeNTOutput ($head, $seq, $out1, $freqStr, $configStr, @OrderArray);
			
			if ($revCom eq "yes")  #running the reverse complement if relevant
			{
				$rev = RevSeq ($chrUse, $head, $seq); #printing the legend inside the sub!
				$head .= "_revcom";
				ElemeNTOutput ($head, $rev, $out1, $freqStr, $configStr, @OrderArray);
				$rev = "";
			}
			$head = "";
			$name = "";
		}

		$seq = "";
		$head = $_;

		#extract relevant input fields to print them out 
			if ($chrUse eq "yes")
			{
				#>hg19_ct_UserTrack_3545_HIST2H2BA range=chr1:120905983-120906083 5'pad=0 3'pad=0 strand=+ repeatMasking=none
				# name = $1;chr = $2;start = $3;end = $4;strand = $5;

				/>(.*)\s.*(chr.+):(\d+)-(\d+).*strand=([-+])/;
				PrintOut::input($2, $3, $4, $5, $1);
			}
			elsif ($chrUse eq "no")
			{
				/>(\S+)/; #printing just the name, no chrs- until the first space
				PrintOut::inputNoChr($1);
			}

	}
	else
	{	
		$_ =~ s/ //g;
		$seq .= $_;
		$seq =~ s!\\[rn]?!!g;  #GC_04_2023 to remove \r and \n
	}
}
	if (length $seq != 0) 
		{ 
			ElemeNTOutput ($head, $seq, $out1, $freqStr, $configStr, @OrderArray); 
			
			if ($revCom eq "yes")  #running the reverse complement if relevant
			{
				$rev = RevSeq ($chrUse, $head, $seq); #printing the legend inside the sub!
				$head .= "_revcom";
				ElemeNTOutput ($head, $rev, $out1, $freqStr, $configStr, @OrderArray);
				$rev = "";
			}
			$head = "";
			$name = "";
		}


sub RevSeq
{
	my ($chrUse, $head, $seq) = @_;
	
	if ($chrUse eq "yes") 
	{
		$head =~ />(.*)\s.*(chr.+):(\d+)-(\d+).*strand=([-+])/;
		$name = $1;
		$name .= "_revcom";
		PrintOut::input($2, $3, $4, $5, $name);
	}
	elsif ($chrUse eq "no")
	{
		$head =~ />(\S+)/; #printing just the name, no chrs- until the first space
		$name = $1;
		$name .= "_revcom";
		PrintOut::inputNoChr($name);
	}			
	
	#print "seq: $seq\n";
	$rev = reverse $seq;
	#print "rev: $rev\n";
	$rev =~ tr/ATGCatgc/TACGtacg/;
	#print "revcom: $rev\n";
	return $rev;
}


sub ElemeNTOutput
{
	
	my ($head, $seq, $out1, $freqStr, $configStr, @OrderArray) = @_;
	$seq =~ s/(\r)|(\n)//g;  #GC_04_2023
	if (length $seq > 1000) 
	{
		PrintOut::PrintSequenceTooLong($head, $out1);
		die;
	}
	
	if ($seq =~ /([^ACGT])/i) #check for non-ACGT characters
	{ 
		my $char = $1;
		PrintOut::ProblemChar($out1, $head, $char);
		die;
	}
	
	
	@ElemeNTSearch::hInrRes = ();
	@ElemeNTSearch::dInrRes = ();
	@ElemeNTSearch::TATARes = ();
	@ElemeNTSearch::BREuRes = ();
	@ElemeNTSearch::BREdRes = ();
	@ElemeNTSearch::hTCTRes = ();
	@ElemeNTSearch::dTCTRes = ();
	@ElemeNTSearch::XCPE1Res = ();
	@ElemeNTSearch::XCPE2Res = ();
	@ElemeNTSearch::PBRes = ();
	@ElemeNTSearch::GAGARes = ();
	#@ElemeNTSearch::IIARERes = ();
	#@ElemeNTSearch::DPRRes = ();
	@ElemeNTSearch::BBCABWRes = ();
	@ElemeNTSearch::Motif1Res = ();

	foreach (@OrderArray)	
	{
		my $elem;
		if ($_ eq 'hInr') {$elem = $ElemeNTSearch::hInr; }
		elsif ($_ eq 'dInr') {$elem = $ElemeNTSearch::dInr; }
		elsif ($_ eq 'TATA') {$elem = $ElemeNTSearch::TATA; }
		elsif ($_ eq 'hTCT') {$elem = $ElemeNTSearch::hTCT; }
		elsif ($_ eq 'dTCT') {$elem = $ElemeNTSearch::dTCT; }
		elsif ($_ eq 'BREu') {$elem = $ElemeNTSearch::BREu; }
		elsif ($_ eq 'BREd') {$elem = $ElemeNTSearch::BREd; }
		elsif ($_ eq 'XCPE1') {$elem = $ElemeNTSearch::XCPE1; }
		elsif ($_ eq 'XCPE2') {$elem = $ElemeNTSearch::XCPE2; }
		elsif ($_ eq 'PB') {$elem = $ElemeNTSearch::PB; } # 5.4.20
		elsif ($_ eq 'GAGA') {$elem = $ElemeNTSearch::GAGA; }
		#elsif ($_ eq 'IIARE') {$elem = $ElemeNTSearch::IIARE; }
		elsif ($_ eq 'DPR') {$elem = $ElemeNTSearch::DPR; }
		elsif ($_ eq 'BBCABW') {$elem = $ElemeNTSearch::BBCABW; }
		elsif ($_ eq 'Motif1') {$elem = $ElemeNTSearch::Motif1; }

		ElemeNTSearch::ElemeNT ($head, $out1, $seq, $elem, $freqStr, $configStr, $gc, 0, 0, 0);
	}
	ElemeNTSearch::GetDependencies($configStr, $seq, $head, $gc);
	ElemeNTSearch::PrintTable($seq, @OrderArray);
	$num ++; #need to acount for die- not represented here
	PrintOut::logSeq ($head);
}

PrintOut::logNum ($num);
PrintOut::Cite; #also print done statement to the terminal


close $fh;
close $out1;
close $log; 
	