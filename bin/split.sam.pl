#!/usr/bin/perl

if($#ARGV != 3)
{
	die("usage: ./split.sam.pl <sam-file> <seq.id.map> <gene.map> <output-dir>\n");
}

my %seqidmap;
open(FILE1, '<', $ARGV[1]) or die("open expression file $ARGV[1] error\n");
while(<FILE1>)
{
	chomp;
	my @line = split('\t');
	$seqidmap{$line[0]} = $line[1];
}
close(FILE1);

print("finish reading seqidmap\n");

# load gene map
my %genemap;
open(FILE1, '<', $ARGV[2]) or die("open expression file $ARGV[2] error\n");
while(<FILE1>)
{
	chomp;
	my @line = split(' ');
	$genemap{$line[0]} = $line[1];
}
close(FILE1);

print("finish reading genemap\n");

# split
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

my %data;
my $pret = "SHAOMINGFU";
my $trunk;

`mkdir -p $ARGV[3]`;

while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	my $tid = $seqidmap{$line[0]};

	if($tid eq $pret)
	{
		$trunk = $trunk . "$_\n";
	}
	else
	{
		if(defined($genemap{$pret}))
		{
			my $gid = $genemap{$pret};
			my $file = $ARGV[3] . "/" . $gid . ".sam";
			open(FILE, '>>', $file) or die("open file $file error\n");
			print FILE $trunk;
			close(FILE);
		}
		$pret = $tid;
		$trunk = "$_\n";
	}
}
