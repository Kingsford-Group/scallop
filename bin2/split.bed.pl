#!/usr/bin/perl

if($#ARGV != 2)
{
	die("usage: ./split.bed.pl <bed-file> <id.map> <output-dir>\n");
}

# load id map
my %idmap;
open(FILE1, '<', $ARGV[1]) or die("open expression file $ARGV[1] error\n");

while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	$idmap{$line[0]} = $line[1];
}

# split
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

my %data;
my $pret = "SHAOMINGFU";
my $trunk;

`mkdir -p $ARGV[2]`;

while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	my @s9 = split(':', $line[3]);
	my $tid = $s9[2];

	if($tid eq $pret)
	{
		$trunk = $trunk . "$_\n";
	}
	else
	{
		if(defined($idmap{$pret}))
		{
			my $gid = $idmap{$pret};
			my $file = $ARGV[2] . "/" . $gid;
			open(FILE, '>>', $file) or die("open file $file error\n");
			print FILE $trunk;
			close(FILE);
		}
		$pret = $tid;
		$trunk = "$_\n";
	}
}
