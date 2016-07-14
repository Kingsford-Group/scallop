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
		$p0[$k] = $p0[$k] + $x[1];
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
	
	next; #TODO

	if($k >= 9)
	{
		my $i = 9;
		$t0[$i] = $t0[$i] + $x[1];
		$t1[$i] = $t1[$i] + $x[5];
		if($x[9] eq "EQUAL")
		{
			$n0[$i] = $n0[$i] + 1;
		}
		if($x[9] eq "LESS")
		{
			$n1[$i] = $n1[$i] + 1;
		}
		if($x[9] eq "GREATER")
		{
			$n2[$i] = $n2[$i] + 1;
		}
	}
}

for (my $i = 1; $i <= 9; $i++)
{
	my $j = $i + 1;
	print("$j $n0[$i] $n1[$i] $n2[$i] $t0[$i] $t1[$i]\n");
}
