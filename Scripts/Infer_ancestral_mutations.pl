#!~/.conda/envs/bioperl/bin/perl

## This script is used to:
## Infer the ancestral mutations.

use warnings;
use strict;
use Bio::TreeIO;
use List::Util qw/max min/;
#mannually set tb.ancestor or tb as the ancestral

die "usage:perl $0 <iqtree_tree_file> <locus_file> 
<ancestor_seq.file> <raw.fasta> <output_db_file> <output_homoplasy_file>\n" if @ARGV==0;

my $treefile = $ARGV[0];
my $treeio = new Bio::TreeIO(-format=>'newick',-file=>$treefile);
my (%ancestor_descendant, %strain_node, %node_sequence, %hash_all, %hash_descendant, $total_node);
if (my $tree = $treeio->next_tree) {
    my @nodes = $tree->get_nodes;
    for (@nodes){
        my $des = $_;
        my $anc = $des->ancestor;
        my $outdes = $des->id;
        if (defined $anc){
            my $outanc = $anc->id;
            if ($outanc =~ /.+/ and $outdes =~ /.+/){
                if ($outdes !~ /tb/){
                    push @{$ancestor_descendant{$outanc}}, $outdes;
                    $hash_descendant{$outdes} = $outanc;
                    $hash_all{$outanc} = 1;
                    $hash_all{$outdes} = 2;
                }
            }
        }
    }
    $total_node = @nodes;
}

"infer_ancestral_mutations.pl" [readonly] 238L, 7367B^[[I            1,1           Top