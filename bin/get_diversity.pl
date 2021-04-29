#!/usr/bin/env perl

use Getopt::Long;
use warnings;
use strict;
use FindBin;

my ($shared, $cpus, $label);

$cpus = 8;
$label = "asv";
GetOptions(
    "s=s" => \$shared,
    "t=i" =>\$cpus,
    "l=s" =>\$label
    );

($shared) ||
    die "usage: $0 OPTIONS\n".
    "where options are:\n".
    "-s <shared file>\n".
    "-t <number of threads>\n".
    "-l <label in shared file: asv|0.03.., default: asv >\n".
    "Read in mothur format shared file to generate alpha, beta diversity\n";


my $bin = "$FindBin::RealBin";
my $projDir = $bin;
$projDir =~ s/\/bin//;
my $mothur = "$projDir/programs/mothur/mothur";


my ($prefix) = $shared =~ /(.*)\.[^.]+$/;
my $prefix_sub = "$prefix.subsample";

#my $sub_sharedFile = "$prefix.$label.subsample.shared";

#rarefaction.single: generate rarefaction curve, rabund files;
#summary.single: get calculator value for each line in the list data and for all possible comparisons between the different groups in the group file
my $cmd = "$mothur \"#rarefaction.single(shared=$shared);summary.single(shared=$shared);";


# summary.shred: for all possible comparisons between the different groups in the group file, produce *.summary file
#Run beta diversity commands
$cmd .= "summary.shared(shared=current, processors=$cpus);get.relabund(shared=current);dist.shared(shared=current, calc=jclass, processors=$cpus);pcoa(phylip=current); nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current,calc=jest, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=braycurtis, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T);dist.shared(shared=current, calc=thetayc, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);\"";
#$cmd .= "otu.association(shared=current, method=spearman, cutoff=0.05);rename.file(input=$prefix.$label.spearman.otu.corr, new=$prefix.$label.spearman.corr);otu.association(shared=current, method=pearson, cutoff=0.05);rename.file(input=$prefix.$label.pearson.otu.corr, new=$prefix.$label.pearson.corr);otu.association(shared=current, method=kendall, cutoff=0.05);rename.file(input=$prefix.$label.kendall.otu.corr, new=$prefix.$label.kendall.corr);\"";
print STDERR "$cmd\n";
system $cmd;


#amova: analyis of molecular variance, to see whether the variation in each sample differs from the variation of the pooled samples
#AMOVA determines whether the genetic diversity within two or more communities is greater than their pooled genetic diversity
#Those P-values less than 0.05 were considered significant.
#Analysis of molecular variance (amova) was applied to test whether the spatial separation of the defined groups visualized in the NMDS plot was statistically significant.
#get diversity for sub sampling
$cmd = "$mothur \"#sub.sample(shared=$shared);rename.file(shared=current, prefix=$prefix.subsample);rarefaction.single(shared=current, processors=$cpus);summary.single(shared=current);summary.shared(shared=current, processors=$cpus);get.relabund(shared=current);dist.shared(shared=current, calc=jclass, processors=$cpus);pcoa(phylip=current); nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current,calc=jest, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=braycurtis, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);dist.shared(shared=current, calc=thetayc, processors=$cpus);pcoa(phylip=current);nmds(phylip=current);amova(phylip=current, design=design.txt);homova(phylip=current, design=design.txt);tree.shared(phylip=current);parsimony(tree=current, group=design.txt, groups=all, processors=$cpus);unifrac.weighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);unifrac.unweighted(tree=current, group=design.txt, groups=all, random=T, processors=$cpus);\"";

#$cmd .= "otu.association(shared=current, method=spearman, cutoff=0.05);rename.file(input=$prefix_sub.$label.spearman.otu.corr, new=$prefix_sub.$label.spearman.corr);otu.association(shared=current, method=pearson, cutoff=0.05);rename.file(input=$prefix_sub.$label.pearson.otu.corr, new=$prefix_sub.$label.pearson.corr);otu.association(shared=current, method=kendall, cutoff=0.05);rename.file(input=$prefix_sub.$label.kendall.otu.corr, new=$prefix_sub.$label.kendall.corr);\"";

print STDERR "$cmd\n";
system $cmd;


