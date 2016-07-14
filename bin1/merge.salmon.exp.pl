#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./merge.salmon.exp.pl <gtf-file> <salmon-expression-file>\n");
}

# load expression file

my %expr;

open(FILE1, '<', $ARGV[1]) or die("open expression file $ARGV[1] error\n");

while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	if($line[1] eq 0)
	{
		next;
	}
	$expr{$line[0]} = $line[1];
}

# read and write gtf file
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	my @s9 = split(' ', $line[8]);
	my $id = substr($s9[3], 1, -2);
	
	if(defined($expr{$id}))
	{
		my $k = index($_, "expression");
		if ($k == -1)
		{
			print("$_ expression \"$expr{$id}\";\n");
		}
		else
		{
			my $ss = substr($_, 0, $k - 1);
			print("$ss expression \"$expr{$id}\";\n");
		}
	}
}
