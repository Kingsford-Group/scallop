#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./add.sam.header.pl <sam-file> <genome-file>\n");
}

my %genome;
open(FILE1, '<', $ARGV[1]) or die("open expression file $ARGV[1] error\n");
while(<FILE1>)
{
	chomp;
	my @line = split('\t');
	$genome{$line[0]} = $line[1];
}
close(FILE1);

# split
my $chr;
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");
while(<FILE2>)
{
	chomp;
	my @line = split('\t');
	$chr = $line[2];
	break;
}
close(FILE2);

print("\@SQ\tSN:$chr\tLN:$genome{$chr}\n");
open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");
while(<FILE2>)
{
	print("$_");
}
close(FILE2);
