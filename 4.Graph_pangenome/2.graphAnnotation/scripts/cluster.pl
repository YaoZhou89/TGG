#!/usr/bin/perl -w
#================================================
#
#         Author: baozhigui
#          Email: zhigui.bao@gmail.com
#         Create: 2021-03-17 23:47:05
#    Description: -
#
#================================================
use strict;
die "\nusage:perl $0 <sort.bed> \n\n"unless (@ARGV==1);
my ($bed)=@ARGV;


my %hash;

my @out;
open I,"$bed" or die $!;
while(<I>)
{
    chomp;
    next if (/^0/);
    my ($chr,$start,$end,$n)=(split /\t/,$_)[0,1,2,3]; 
    if ($n =~ /TomatoPan/)
    {

    }else{
        my ($gc)=$n=~/.*_(\d+)T.*/;
        $gc=~s#^0##g;
        next if ($gc ne $chr);
    }
    next if ($end -$start > 20000);
    #print "$_\n";
    my $key=$chr."\t".$start."\t".$end;
    if (exists $hash{$key})
    {
        next;
    }else{
        push @out,$key;
        $hash{$key}=$_;
    }
}
close I;

foreach my $k(@out)
{
    print "$hash{$k}\n";
}
