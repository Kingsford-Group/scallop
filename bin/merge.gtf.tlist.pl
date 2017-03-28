#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./merge.gtf.tlist.pl <gtf-file> <tlist-file>\n");
}

# load expression file

my %expr;

open(FILE1, '<', $ARGV[1]) or die("open expression file $ARGV[1] error\n");

while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	$expr{$line[0]} = $line[3];
}

# read and write gtf file
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	my @s9 = split(' ', $line[8]);
	my $id = substr($s9[1], 1, -2);
	
	if(defined($expr{$id}))
	{
		print("$_ expression \"$expr{$id}\";\n");
	}
}
