#!/usr/bin/perl

if($#ARGV != 0)
{
	die("usage: ./compare.gtf.pl <gtf-file>\n");
}

open(FILE, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");
while(<FILE>)
{
	chomp;
	my @line = split('\t');
		
	$f = 0;
	for($k = 1; $k <= 22; $k++)
	{
		if($line[0] eq "chr$k")
		{
			$f = 1;
			last;
		}
	}

	if($line[0] eq "chrX")
	{
		$f = 1;
	}

	if($line[0] eq "chrY")
	{
		$f = 1;
	}

	if($f eq 1)
	{
		printf("$_\n");
	}
}
