#!/usr/bin/perl -w
#================================================
#
#         Author: baozhigui
#          Email: zhigui.bao@gmail.com
#         Create: 2021-03-17 15:33:46
#    Description: -
#
#================================================
use strict;
die "\nusage:perl $0 <gmap_gff3> <coverage> <identity>\n\n"unless (@ARGV==3);
my ($g,$c,$i)=@ARGV;

my @name;
my %hash;
open I,"$g" or die $!;
while(<I>)
{
    chomp;
    next if (/^#/);
    my @ll = (split /\t/,$_);
    next if ($ll[2] eq "gene");
    if ($ll[2] eq "mRNA")
    {
        my ($id,$cov,$iden)=$ll[8]=~/^ID=(.*?)\.mrna1;Name.*coverage=(.*?);identity=(.*?);/;
        
        if ($cov > $c and $iden > $i)
        {
            $hash{$id}="$_";
            push @name,$id;
        }
        #print "$id\t$cov\t$iden\n";
    }else{
        my ($pid) = $ll[8]=~/^ID=(.*?)\.mrna1.*;/;
        $hash{$pid}.="\n$_";
    }
}
close I;

foreach my $k (@name)
{
    print "$hash{$k}\n";
}
