#!/usr/bin/perl

if($#ARGV != 0)
{
	die("usage: ./summary.result.pl <all.cmp>\n");
}

# load expression file

my @n0;	# EQUAL instances
my @n1; # LESS instances
my @n2; # GREATER instances
my @t0; # total predicted transcripts
my @t1; # total correct transcripts

open(FILE, '<', $ARGV[0]) or die("open expression file $ARGV[0] error\n");

for (my $i = 0; $i <= 12; $i++)
{
	$n0[$i] = $n1[$i] = $n2[$i] = $t0[$i] = $t1[$i] = 0;
}

while(<FILE>)
{
	chomp;
	my @x = split(' ');
	my $k = $x[3] - 1;

	if($k <= 9)
	{
		$t0[$k] = $t0[$k] + $x[1];
		$t1[$k] = $t1[$k] + $x[5];
		if($x[9] eq "EQUAL")
		{
			$n0[$k] = $n0[$k] + 1;
		}
		if($x[9] eq "LESS")
		{
			$n1[$k] = $n1[$k] + 1;
		}
		if($x[9] eq "GREATER")
		{
			$n2[$k] = $n2[$k] + 1;
		}
	}
#if( ($k >= 10) and ($k <= 11) )
	if($k >= 10)
	{
		$t0[10] = $t0[10] + $x[1];
		$t1[10] = $t1[10] + $x[5];
		if($x[9] eq "EQUAL")
		{
			$n0[10] = $n0[10] + 1;
		}
		if($x[9] eq "LESS")
		{
			$n1[10] = $n1[10] + 1;
		}
		if($x[9] eq "GREATER")
		{
			$n2[10] = $n2[10] + 1;
		}
	}
	if( ($k >= 12) and ($k <= 14) )
	{
		$t0[11] = $t0[11] + $x[1];
		$t1[11] = $t1[11] + $x[5];
		if($x[9] eq "EQUAL")
		{
			$n0[11] = $n0[11] + 1;
		}
		if($x[9] eq "LESS")
		{
			$n1[11] = $n1[11] + 1;
		}
		if($x[9] eq "GREATER")
		{
			$n2[11] = $n2[11] + 1;
		}
	}
	if($k >= 15)
	{
		$t0[12] = $t0[12] + $x[1];
		$t1[12] = $t1[12] + $x[5];
		if($x[9] eq "EQUAL")
		{
			$n0[12] = $n0[12] + 1;
		}
		if($x[9] eq "LESS")
		{
			$n1[12] = $n1[12] + 1;
		}
		if($x[9] eq "GREATER")
		{
			$n2[12] = $n2[12] + 1;
		}
	}
}

for (my $i = 1; $i <= 10; $i++)
{
	my $j = $i + 1;
	print("$j: $n0[$i] $n1[$i] $n2[$i] $t0[$i] $t1[$i]\n");
}
