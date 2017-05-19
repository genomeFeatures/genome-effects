#!/usr/bin/perl

#exon count vs Transcript count by organism version and gene prediction id
$dbname="graber_transcriptdb"; $host="harlequin"; $user="pup";$pass="puppass";
$dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;mysql_local_infile=1",$user, $pass);

if(!$dbh){print  "Could not connect to database :$!\n"; exit;} 
 
$query="select gene_prediction_id, transcript_name,count(distinct exon_id) 
        from transcript_exon where organism_version_id=? group by  gene_prediction_id, transcript_name";
$qh_getExFreq=$dbh->prepare($query);

$query="select gene_prediction_id, transcript_id,length(tx_end-tx_start) 
        from transcript t, transcript_by_annotation  ta  
        where t.organism_version_id=? and t.organism_version_id=ta.organism_version_id ";
$qh_getTxLen=$dbh->prepare($query);



