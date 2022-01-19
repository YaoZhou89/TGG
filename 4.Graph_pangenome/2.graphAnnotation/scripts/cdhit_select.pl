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
die "\nusage:perl $0 <cd-hit.clstr> <order.txt>\n\n"unless (@ARGV==2);
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
$/=">Cluster ";
open I,"$bed" or die $!;
<I>;
while(<I>)
{
    chomp;
    #print "$_\n";
    my ($key,$seq)=(split /\n/,$_,2)[0,1];
    #print "$key\n";
    my @tmp = (split /\n/,$seq);
    foreach my $k(@tmp)
    {
        my $g=(split /,/,$k)[1];
        #print "$name\n";
        $g=~s# >##g;
        my $name=(split /\.\.\./,$g)[0];

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
                $gene{$key}=$name;
                $hash{$key}=$ord;
            }
        }else{
            $gene{$key}=$name;
            $hash{$key}=$ord;
        }

    }

}
close I;

#print Dumper \%hash;
foreach my $k (keys %gene)
{
    print "$gene{$k}\t$hash{$k}\n";
}
