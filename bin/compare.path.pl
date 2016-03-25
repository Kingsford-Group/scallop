#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./compare.count.pl <file1> <file2>\n");
}

my %m1;;;;
open(FILE1, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");
while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	$m1{$line[0]} = $line[1];
}

my %m2;
open(FILE2, '<', $ARGV[1]) or die("open gtf file $ARGV[1] error\n");
while(<FILE2>)
{
	chomp;
	my @line = split(' ');
	$m2{$line[0]} = $line[1];
}

foreach my $s (keys %m2)
{
	if(defined($m1{$s}))
	{
		if($m1{$s} < $m2{$s})
		{
			print("$s $m1{$s} $m2{$s} BETTER\n");
		}
		elsif($m1{$s} > $m2{$s})
		{
			print("$s $m1{$s} $m2{$s} WORSE\n");
		}
		else
		{
			print("$s $m1{$s} $m2{$s} EQUAL\n");
		}
	}
}
