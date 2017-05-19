#!/usr/bin/perl
#Variations in the genome
#One of the most important uses for the variation database is to be able to get all variations in a certain region in the genome. Below it is a simple commented perl script to illustrate how to get all variations in chromosome 25 in zebrafish. 
#
#
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

my $slice_adaptor = $registry->get_adaptor('danio_rerio', 'core', 'slice'); #get the database adaptor for Slice objects
my $slice = $slice_adaptor->fetch_by_region('chromosome',25); #get chromosome 25 in zebrafish

my $vf_adaptor = $registry->get_adaptor('danio_rerio', 'variation', 'variationfeature'); #get adaptor to VariationFeature object
my $vfs = $vf_adaptor->fetch_all_by_Slice($slice); #return ALL variations defined in $slice

foreach my $vf (@{$vfs}){
  print "Variation: ", $vf->variation_name, " with alleles ", $vf->allele_string, 
        " in chromosome ", $slice->seq_region_name, " and position ", $vf->start,"-",$vf->end,"\n";
}
