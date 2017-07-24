use graber_transcriptdb;

select distinct gene_prediction_name,sample_name,ta.transcript_name,ga.gene_name 
from rnaseq_sample s, rnaseq_transcript rt, gene_prediction p, transcript_by_annotation ta, gene_by_annotation ga 
where s.sample_id=rt.sample_id 
and rt.transcript_name=ta.transcript_name 
and rt.gene_prediction_id=ta.gene_prediction_id 
and rt.organism_version_id=ta.organism_version_id 
and ta.transcript_id=ga.transcript_id 
and ta.gene_prediction_id=ga.gene_prediction_id 
and ga.gene_prediction_id=p.gene_prediction_id 
and p.gene_prediction_id=rt.gene_prediction_id;
