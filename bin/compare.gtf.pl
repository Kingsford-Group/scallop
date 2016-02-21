#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./compare.gtf.pl <gtf-file1> <gtf-file2>\n");
}

my %m1;
open(FILE1, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");
while(<FILE1>)
{
	chomp;
	my @line = split('\t');

	if($line[2] ne "transcript")
	{
		next;
	}

	my @s9 = split(' ', $line[8]);

	my $k;
	for($k = 0; $k < length(@s9); $k++)
	{
		if($s9[$k] eq "gene_id")
		{
			last;
		}
	}

	my $gid = substr($s9[$k + 1], 1, -2);
	$m1{$gid} = $m1{$gid} + 1;
}

my %m2;
open(FILE2, '<', $ARGV[1]) or die("open gtf file $ARGV[1] error\n");
while(<FILE2>)
{
	chomp;
	my @line = split('\t');

	if($line[2] ne "transcript")
	{
		next;
	}

	my @s9 = split(' ', $line[8]);

	my $k;
	for($k = 0; $k < length(@s9); $k++)
	{
		if($s9[$k] eq "gene_id")
		{
			last;
		}
	}

	my $gid = substr($s9[$k + 1], 1, -2);
	$m2{$gid} = $m2{$gid} + 1;
}

foreach my $s (keys %m1)
{
	if(defined($m2{$s}))
	{
		if($m1{$s} > $m2{$s})
		{
			print("$s $m1{$s} $m2{$s} BETTER\n");
		}
		elsif($m1{$s} < $m2{$s})
		{
			print("$s $m1{$s} $m2{$s} WORSE\n");
		}
		else
		{
			print("$s $m1{$s} $m2{$s} EQUAL\n");
		}
	}
	else
	{
		print("$s $m1{$s} 0 FALSE\n");
	}
}

foreach my $s (keys %m2)
{
	if(defined($m1{$s}))
	{
		next;
	}
	print("$s 0 $m2{$s} FALSE\n");
}
