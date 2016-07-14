#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./map.salmon.exp.pl <ensembl-name-file> <salmon-sparse-quant-file>\n");
}

# load expression file

open(FILE1, '<', $ARGV[0]) or die("open ensembl-name file $ARGV[0] error\n");

my @list;
while(<FILE1>)
{
	chomp;
	if ($line[0] eq '#')
	{
		next;
	}

	my @line = split('\|');
	if($#line < 1 or $line[1] eq "" or $line[1] !~ /ENST/)
	{
		$list[$cnt] = "ABCDEFG";
	}
	else
	{
		$list[$cnt] = $line[1];
	}
	$cnt++;
}

open(FILE2, '<', $ARGV[1]) or die("open salmon-sparse-quant-file $ARGV[1] error\n");

while(<FILE2>)
{
	chomp;
	if ($line[0] eq '#')
	{
		next;
	}

	my @line = split(',');
	if($#line < 1 or $line[1] eq "")
	{
		next;
	}

	my $k = $line[0];
	if($k > $#list or $list[$k] eq "ABCDEFG")
	{
		next;
	}

	my $exp = int($line[1] * 10);

	if($exp <= 9)
	{
		next;
	}

	print("$list[$k] $exp\n");
}
