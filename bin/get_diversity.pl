#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;
use FindBin;

my ($list, $groups, $cpus, $mothur);
$cpus = 1;
GetOptions(
    "l=s" => \$list,
    "g=s" => \$groups,
    "t=i" =>\$cpus,
    "p=s" =>\$mothur
    );

($list && $groups) ||
    die "usage: $0 OPTIONS\n".
    "where options are:\n".
    "-l <list file>\n".
    "-g <group file>\n".
    "-t <number of threads>\n".
    "-p <mothur program path>\n".
    "Read in mothur list and group file to generate alpha, beta diversity\n";


my $bin = "$FindBin::RealBin";
#my $cmd = "$bin/get_diversity.pl -l $analysisname.list -g $analysisname.groups -t $cpus ";

my ($prefix) = $list =~ /(.*)\.[^.]+$/;
my $sharedFile = "$prefix.shared";
my $sub_sharedFile = "$prefix.0.03.subsample.shared";
my $cmd = "$mothur \"#make.shared(list=$list, group=$groups);get.otulist(list=current);";

#rarefaction.single: generate rarefaction curve, rabund files;
#summary.single: get calculator value for each line in the list data and for all possible comparisons between the different groups in the group file
$cmd .= "rarefaction.single(shared=current);summary.single(shared=current);";


# summary.shred: for all possible comparisons between the different groups in the group file, produce *.summary file
#Run beta diversity commands
$cmd .= "summary.shared(shared=current, processors=$cpus);get.relabund(shared=current);dist.shared(shared=current, calc=jclass, processors=$cpus);pcoa(phylip=current); nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current,calc=jest, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=braycurtis, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T);dist.shared(shared=current, calc=thetayc, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);\"";

print STDERR "$cmd\n";
system $cmd;

$cmd = "$bin/filter_shared.pl -s $sharedFile -mp 0.1 -mns 2 > network.shared;";
$cmd .= "$mothur \"#otu.association(shared=network.shared, method=spearman, cutoff=0.05);otu.association(shared=network.shared, method=pearson,cutoff=0.05);otu.association(shared=network.shared, method=kendall, cutoff=0.05);\"";

print STDERR "$cmd\n";
system $cmd;



#amova: analyis of molecular variance, to see whether the variation in each sample differs from the variation of the pooled samples
#AMOVA determines whether the genetic diversity within two or more communities is greater than their pooled genetic diversity
#Those P-values less than 0.05 were considered significant.
#Analysis of molecular variance (amova) was applied to test whether the spatial separation of the defined groups visualized in the NMDS plot was statistically significant.
#get diversity for sub sampling
$cmd = "$mothur \"#sub.sample(list=$list,group=$groups, persample=t);make.shared(list=current, group=current);rarefaction.single(shared=current, processors=$cpus);summary.single(shared=current);summary.shared(shared=current, processors=$cpus);get.relabund(shared=current);dist.shared(shared=current, calc=jclass, processors=$cpus);pcoa(phylip=current); nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current,calc=jest, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=braycurtis, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=thetayc, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);\"";




#otu.association(shared=current, method=spearman, cutoff=0.05);otu.association(shared=current, method=pearson, cutoff=0.05);otu.association(shared=current, method=kendall, cutoff=0.05);\"";
print STDERR "$cmd\n";
system $cmd;

$cmd = "$bin/filter_shared.pl -s $sub_sharedFile -mp 0.1 -mns 2 > network.0.03.subsample.shared;";
$cmd .= "$mothur \"#otu.association(shared=network.0.03.subsample.shared, method=spearman, cutoff=0.05);otu.association(shared=network.0.03.subsample.shared, method=pearson,cutoff=0.05);otu.association(shared=network.0.03.subsample.shared, method=kendall, cutoff=0.05);\"";

print STDERR "$cmd\n";
system $cmd;

