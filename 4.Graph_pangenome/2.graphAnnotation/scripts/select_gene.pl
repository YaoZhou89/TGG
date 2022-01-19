#!/usr/bin/perl -w
#================================================
#
#         Author: baozhigui
#          Email: zhigui.bao@gmail.com
#         Create: 2021-03-18 01:07:15
#    Description: -
#
#================================================
use strict;
use Data::Dumper;
die "\nusage:perl $0 <bed> <order.txt>\n\n"unless (@ARGV==2);
my ($bed,$o)=@ARGV;



my %order;
my $count=0;
open I,"$o" or die $!;
while(<I>)
{
    chomp;
    $count+=1;
    $order{$_}=$count;   
}
close I;




my %hash;
my %gene;
open I,"$bed" or die $!;
while(<I>)
{
    chomp;
    my ($chr,$start,$end,$gs,$ge,$name)=(split /\t/,$_)[0,1,2,4,5,6];
    my $key = $chr."\t".$start."\t".$end;
    my $or;

    if ($name=~/Solyc_cer/)
    {
        ($or)=$name=~/Solyc_cer_(.*?)_/;
    }elsif ($name =~ /Pan/)
    {
        $or = "TomatoPan"
    }else{
        ($or)=$name=~/^S.*_(.*?)_/;
    }

    my $ord = $order{$or};
    if (exists $hash{$key})
    {
        if ($ord < $hash{$key})
        {
            $gene{$key}=$_;
            $hash{$key}=$ord;
        }
    }else{
        $gene{$key}=$_;
        $hash{$key}=$ord;
    }
}
close I;

#print Dumper \%hash;
foreach my $k (keys %gene)
{
    print "$gene{$k}\n";
}
