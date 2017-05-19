<?php
 $server="http://demon.jax.org";
 $genome_list_file="$server/genomes/service_organisms_w_chromosome_data.txt";
 $dir=getcwd();
 $service_base="$server/transcriptdb/web-services/";
 $usage= "
      <div style='width:95%;margin:auto; font-size:0.95em;font-family:\"Helvetica Neue Light\", HelveticaNeue-Light, \"elvetica Neue\", Helvetica, Arial, sans-serif'>
        <style>li{padding:0.4em;}h1{color:#800;}#note{border:dashed 1px #eee;background:#efefef;
            width:600px;margin:auto;padding:0.8em;}p,li{color:#333333;}</style>
        <h1>Graber Transcript Database Web Services</h1>
        <p>
          This web service generates transcript annotations that overlap with a given genomic range 
          or genomic term. A term can be a gene or transcript symbol,or id .
         The result is either in json or xml format(default xml).
         By default, the web service returns the isoform view of the gene. If you call the service
         with the isoform flag set to 0, the service returns the agregate view of the gene.
       </p>
       <div id='note'>
       <b>Note:</b>
       Features are stored in our database in zero-based cordinates standard for feature starts
       and 1-based cordinate standard for feature ends. We display all feature coordinates in 1-base.
       </div>
       <h2>How To Use It:</h2>
         <ol>
          <li> Display this help document: 
            <a href='$service_base?h=1'>$service_base?h=1</a></li>
        <li> Display the list of current gene annotation sources:
              <a href='$service_base?l=1'>$service_base?l=1</a><br/></li>
        <li> Display the list of current organism versions:
              <a href='$service_base?s=1'>$service_base?s=1</a><br/></li>
         <li> Display the list of organism versions with genome sequence:
              <a href='$service_base?g=1'>$service_base?g=1</a><br/></li>
        <li> Display the list of gene annotation sources for mm9:
              <a href='$service_base?v=mm9&l=1'>$service_base?v=mm9&l=1</a><br/></li>
        
       <li>>>query: <a href='$service_base?v=mm9&t=pax6&a=mm9-ensGene&format=json'>$service_base?t=pax6&amp;a=mm9-ensGene&amp;format=json</a><br/>
          <b>description:</b> This will generate mm9 Ensembl annotations that were provided for pax6
             and return the isoform view of the gene in json format <br/></li>

    
     <li>>>query: <a href='$service_base?format=xml&v=mm10&t=4:120405281-120454515&a=GRCm38-mgiGene'>
                  $service_base?format=xml&v=mm10&t=4:120405281-120454515&a=GRCm38-mgiGene</a></li>
      <b>description:</b> Use MGI annotations to return the isoform view of all the transcripts/exons 
                on chromosome 4 that overlap mm10 genomic region 120405281-120454515 on both + and - strands <br/></li>
      <li>>>query: <a href='$service_base?format=xml&v=mm10&t=4:120405281-120454515&a=GRCm38-mgiGene&isoform=0'>
                  $service_base?format=xml&v=mm10&t=4:120405281-120454515&a=GRCm38-mgiGene&isoform=0</a></li>
      <b>description:</b> Use MGI annotations to return the agregate view of all the transcripts/exons 
                on chromosome 4 that overlap mm10 genomic region 120405281-120454515 <br/></li>
       </ol>
       </ol>
     
      <h2>Arguments List:</h2><ul>
      <li>h=1  displays this help page</li>
      <li>l=1  displays the list of current gene prediction sources found in Graber Transcript Database</li>
      <li>s=1  displays the list of current organism versions found in Graber Transcript Database</li>
     <li>t=term is a query term. <b>gene symbol, transcript accession id,gene accession id,<br/>a genomic range</b>
          (chrom:locStart-locEnd) where:<br/>
           a. chrom is the chromosome name (a string, 1,2,3,X,or chr1,chr2,chrX,...)<br/>
           b. locStart/locEnd are starting and ending positions in bp or in mbp (40.5)
            </li>
     <li>strand=[+/-] strand is the genome strand (+,-) is optional, the default is to generate hits on both (+ and -)</li>
     <li>format=output_format  is the output format, xml or json (default standard xml)</li>
     <li>v=organism_version  is the ucsc organism version (default mm9)</li>
     <li> a=gene_prediction_name  gene prediction name <br/>
       <p>Note: Every organism has the default gene prediction, for queries using
          the t option (search term for example a gene/transcript name,or id),if the query term was not provided by the
          defaullt or the specified gene prediction, we set the gene prediction
          to be what we have in the database for that term
       </p>
     </li>
    <li>isoform=[1/0]; if isoform=1 returns gene result one 
           isoform per row(default), if isoform=0 returns aggregate transcript and exons</li>
    </ul>
   <h2>Result structure and fields:</h2>
    a. annotations: document root<br/>
    b. head : contains the source name, organism version used to generate the annotations<br/>
    c. annotation: a row representation of every item found. The following fields are includes:<br/>
     <p>c.1  chrom: is the chromosome name <br/>
     c.2  chromStart : gemomic location start (relative to the query term)<br/>
     c.3  chromEnd : gemomic location end (relative to the query term)<br/>
     c.4  strand : genomic strand<br/>
     c.5  geneNames: gene list that overlap your query<br/>
     c.6. annotationStart : genomic location start (relative to overlap annotations)<br/>
     c.7  annotationEnd : genomic location end (relative to overlap annotations)<br/>
     c.8  transcriptNames:transcript (list) that overlap your query <br/>
     c.9  exonStarts: list of exon starts that overlap your query<br/>
     c.10 exonEnds : list of exon ends that overlap your query
   
   </div>
  ";
?>
