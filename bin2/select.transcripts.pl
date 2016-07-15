#!/usr/bin/perl

if($#ARGV != 1)
{
	die("usage: ./select.transcripts.pl flux-profile id.map\n");
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

my %count;
my @list;
while(<FILE2>)
{
	chomp;
	my @s = split('\t');
	my $tid = $s[1];

	# filter transcripts
	if($s[5] < 10 or $s[7] <= 0 or $s[9] <= 0 or $s[10] < 0.8)
	{
		next;
	}

	if(defined($idmap{$tid}))
	{
		my $gid = $idmap{$tid};
		if(defined($count{$gid}))
		{
			$count{$gid} = $count{$gid} + 1;
		}
		else
		{
			$count{$gid} = 1;
		}

		my $x = "$gid $tid $s[2] $s[3] $s[5] $s[7] $s[9] $s[10]";
		push(@list, $x);
#print("$x\n");
	}
}

foreach my $k (@list)
{
	my @s = split(' ', $k);
	if($count{$s[0]} >= 3)
	{
#print("$k\n");
		print("$s[1] $s[0]\n");
	}
}
