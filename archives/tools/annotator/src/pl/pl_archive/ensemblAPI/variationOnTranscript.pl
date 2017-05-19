#!/usr/bin/perl
#In the example below, it is explained how to get all variations in a particular human transcript and see what is the effect of that 
#variationin the transcript, including the PolyPhen and SIFT predictions. 
#It is also shown how to retrieve the Sequence Ontology terms for the consequences.
#http://useast.ensembl.org/info/docs/api/variation/variation_tutorial.html#annotation
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);
my  $dba =$registry->get_DBAdaptor("Mouse", "core");
my $stable_id = 'ENSMUST00000107826'; #this is the stable_id of a human transcript
my $transcript_adaptor = $registry->get_adaptor('mus_musculus', 'core', 'transcript'); 
#my $db_connection = $dba->dbc();
q{printf(
        "species/group\t%s/%s\ndatabase\t%s\nhost:port\t%s:%s\n\n",
        $dba->species(),  $dba->group(),
        $db_connection->dbname(), $db_connection->host(),
        $db_connection->port()
    );
};
#exit(0);
my $transcript = $transcript_adaptor->fetch_by_stable_id($stable_id); #get the Transcript object
my $trv_adaptor = $registry->get_adaptor('mus_musculus', 'variation', 'transcriptvariation'); #get the adaptor to get TranscriptVariation objects
# Fetch a variation object
my $va_adaptor = $registry->get_adaptor('mus_musculus', 'variation', 'variationannotation');
my $var_adaptor = $registry->get_adaptor('mouse', 'variation', 'variation');
my $trvs = $trv_adaptor->fetch_all_by_Transcripts([$transcript]); #get ALL effects of Variations in the Transcript
foreach my $tv (@{$trvs}) {
	my $tvas = $tv->get_all_alternate_TranscriptVariationAlleles();
	foreach my $tva(@{$tvas}) {
		my @ensembl_consequences;my @so_consequences;
		my $ocs = $tva->get_all_OverlapConsequences();
		foreach my $oc(@{$ocs}) {
			push @ensembl_consequences, $oc->display_term;
			push @so_consequences, $oc->SO_term;
		}
		my $sift = $tva->sift_prediction;
		my $polyphen = $tva->polyphen_prediction;
		print "Variation ", $tv->variation_feature->variation_name," allele ", $tva->variation_feature_seq,
                       " ,Alleles ", $tva->allele_string," ,Chrom: ",$tva->seq_region_name, ":", $tva->start,"-",$tva->end,
		       " has consequence ", join(",", @ensembl_consequences)," (SO ", join(",", @so_consequences), ")\n";
	} 
}

