#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./split.gtf.pl <gtf-file> <output-dir>\n");
}

# split
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

my $pret = "SHAOMINGFU";
my $trunk;

`mkdir -p $ARGV[1]`;

while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	my @s9 = split(' ', $line[8]);
	my $gid = substr($s9[1], 1, -2);

	if($gid eq $pret)
	{
		$trunk = $trunk . "$_\n";
	}
	else
	{
		if($pret ne "SHAOMINGFU")
		{
			my $file = $ARGV[1] . "/" . $pret . ".gtf";
			open(FILE, '>>', $file) or die("open file $file error\n");
			print FILE $trunk;
			close(FILE);
		}

		$pret = $gid;
		$trunk = "$_\n";
	}
}
