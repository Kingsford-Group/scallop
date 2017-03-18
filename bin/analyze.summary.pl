#!/usr/bin/perl

if($#ARGV != 0)
{
	die("usage: ./analyze.summary.pl <summary-file>\n");
}

open(FILE, '<', $ARGV[0]) or die("open summary file $ARGV[0] error\n");

my $tp1 = 0;
my $tp2 = 0;
my $ss1 = 0;
my $ss2 = 0;
my $ref = 0;

while(<FILE>)
{
	chomp;
	my @x = split(' ');
	$ss1 += $x[1];
	$ss2 += $x[5];
	$ref += $x[2];
	$tp1 += $x[2] * $x[3] / 100;
	$tp2 += $x[6] * $x[7] / 100;
}

my $sn1 = $tp1 / $ref * 100.0;
my $sn2 = $tp2 / $ref * 100.0;
my $sp1 = $tp1 / $ss1 * 100.0;
my $sp2 = $tp2 / $ss2 * 100.0;

my $m = sprintf("reference = %d, prediction = (%d, %d), correct = (%d, %d), sensitivity = (%.2f, %.2f), specificity = (%.2f, %.2f)\n", 
		$ref, $ss1, $ss2, $tp1, $tp2, $sn1, $sn2, $sp1, $sp2);

print("$m");
