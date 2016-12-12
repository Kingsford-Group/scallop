#!/usr/bin/perl

if($#ARGV != 0)
{
	die("usage: ./stat.gtf.pl <gtf-file>\n");
}

open(FILE2, '<', $ARGV[0]) or die("open gtf file $ARGV[0] error\n");

my %gcount;
my %tmap;

while(<FILE2>)
{
	chomp;
	my @line = split('\t');

	if($line[2] ne "transcript")
	{
		next;
	}

	my @s9 = split(' ', $line[8]);
	my $gid = substr($s9[1], 1, -2);
	my $tid = substr($s9[3], 1, -2);

	$tmap{$tid} = $gid;

	if(defined($gcount{$gid}))
	{
		$gcount{$gid} = $gcount{$gid} + 1;
	}
	else
	{
		$gcount{$gid} = 1;
	}
}

for my $k (keys %tmap)
{
	$gid = $tmap{$k};
	$count = $gcount{$gid};
	print("$k $gid $count\n");
}
