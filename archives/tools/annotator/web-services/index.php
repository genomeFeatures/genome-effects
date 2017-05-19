<?php 
/********************** Transcript Database Web Service **************************
 Author: Lucie Ndzana Hutchins
 Date : July 2012
  
  Note: This utility processes user's query to return the appropriate results.
        The user can specify the result format: json or xml - xml is the default
        If called without any argument, the result will be a "help me" page

 Features are stored in our database in zero-based cordinates standard for
 feature starts and 1-based cordinate standard for feature ends.

Arguments List:
    h=1 displays the help page  Example: /web-services/ or Example: /web-services/?h=1
    s=1 returns the list of organism with gene annotations. If "host" is specified, gets the results from
        the specified host. The default our local annotations database
    v=organism_version  is the ucsc organism version (default mm9)
    l=1 returns the list of gene prediction sources for the specified organism version or all listing if no organism version specified
        Example: /web-services/?v=mm9&l=1  will list all gene annotation sources  we have in our db for mm9  
   
    t=term is a query term. A term can be: gene symbol, transcript accession id,gene accession id,
         or a genomic range (chrom:locStart-locEnd) where:
            a. chrom is the chromosome name (a string, 1,2,3,X,or chr1,chr2,chrX,...)
            b. locStart/locEnd are starting and ending positions in bp or in mbp (40.5)
    strand=[+/-] strand is the genome strand (+,-) is optional, the default is to generate hits on both (+ and -)
            strand option is checked only when user queries using a term ( the t option)
    format=output_format is the output format, xml or json (default standard xml)
    a=gene_prediction_name gene prediction source name for example ensGene, mgiGene,refGene
        Example: /web-services/?v=mm9&a=refGene will use refGene annotations for mm9 to run the query
    aggregate=[1/0]; if aggregate=0 returns gene result one isoform per row(default), 
                     else returns aggregate transcript and exons (one gene by row)
     
    host=ucsc specifies the remote database server to use for annotations extraction

  Note: Every organism has the default gene prediction, for queries using the t option 
        (search term for example a gene/transcript name,or id),
        if the query term was not provided by the defaullt or the specified gene prediction, 
        we set the gene prediction to be what we have in the database for that term

    
***********************************************************************************/
 $filename="tx-service.json";$org_vid=0;$pred_id=0;$mb=1000000;
 //set up the result format
 $format = strtolower($_GET['format']) == 'json' ? 'json' : 'xml'; 
 if($format=="xml")$filename="tx-service.xml";
 //set up the remote server
 $host=isset($_GET["host"])?$_GET["host"]:"";
 //get user's selections 
 $term=isset($_GET["t"])?$_GET["t"]:"";
 $strand=isset($_GET['strand'])? $_GET['strand']:"";
 $table_name="";$isoform=isset($_GET['aggregate'])?0:1; 

 $oganism_listing=isset($_GET['s'])? intval($_GET['s']) : 0;
 $annotation_listing=isset($_GET['l'])? intval($_GET['l']) : 0;
 $genome_listing =isset($_GET['g'])? intval($_GET['g']) : 0;
 $chromosome_listing=isset($_GET["chrl"])? intval($_GET["chrl"]) : 0;
 $annotation_load=isset($_GET['annot'])? intval($_GET['annot']) : 0;
 $atype=isset($_GET['atype'])? $_GET['atype']: "";
 $feature_listing=isset($_GET['featureList'])?intval($_GET['featureList']):0;
 //set result file header
 if((empty($term)||($term==""))&&(!$annotation_listing)&&(!$oganism_listing)&&
    (!$genome_listing)&&!($chromosome_listing))$format="text/html";
 if(isset($_GET['h']))$format="text/html";
 if($format=="json" || $format=="xml"){
    header("Content-type: application/octet-stream");
    header("Cache-Control: cache, must-revalidate");
    header("Pragma: public");
    header("Content-type: application/$format");
  }
  require_once("transcript_utils.php");require_once("global.php");global $usage;
  require_once("ucsc_transcript_utils.php");
  $organism_version=isset($_GET['v']) ? filterQuery(trim(strtolower($_GET['v']))) : "";
  $term=filterQuery($term); $strand=filterQuery($strand);$host=filterQuery($host);
  $db_name=($host=="")?"Graber Transcript ":strtoupper($host);
  $org_vid=getOrganismVid($organism_version);
  $query_type=isset($_GET["querytype"])? filterQuery($_GET["querytype"]) : "";
  $gene_prediction=isset($_GET['a'])? filterQuery($_GET['a']): "All";  $title=" $organism_version -"; 

  if($oganism_listing) $title.=" List of organism with gene annotations";
  if($feature_listing)$title.=" List of annotated features for the specified organism assembly version";
  if($annotation_listing) $title.=" List of available gene annotation sources";
  if($genome_listing) $title.=" List of organism with available chromosome assembly data";
  if($annotation_load) $title.=" Generating $gene_prediction annotations file ";
  $head= "Using $db_name database -$title"; $data_desc="Result fiedls:";
  $xml_shebang=' <?xml version="1.0" encoding="ISO-8859-1"\? >';
  $annotations[]=array("head"=>$head);$annotation_list=array();
  //now return result based on query type  
  if(isset($_GET['h'])){if($format=="json")echo json_encode($usage);else echo $usage;}
  elseif($oganism_listing){$db; //the query is to list organisms with gene annotations 
     if(empty($host)){$db=getOrgWithGenome();displayOrganismList($db,$annotation_list);}
     else{ ucsc_displayOrganismList($db,$annotation_list);}
     $annotations=array_merge($annotations,$annotation_list);
     if($format=="json")echo json_encode(array('organisms'=>$annotation_list));
     else{echo'<organisms>'.displayXml($annotations).'</organisms>';}
  }
  elseif($annotation_listing){//the query is to list all gene annotations source for the specified organism assembly version
     if(empty($host))getListGenePrediction($org_vid,$annotation_list);
     else ucsc_getListGenePrediction($organism_version,$annotation_list);
     $annotations=array_merge($annotations,$annotation_list);
     if($format=="json")echo json_encode(array('gene_predictions'=>$annotation_list));
     else{echo'<gene_predictions>'.displayXml($annotations).'</gene_predictions>';}
  }
 elseif($feature_listing){//the query is to list features other than genes for the specified organism assembly version
     ucsc_getListFeatures($organism_version,$annotation_list);
     $annotations=array_merge($annotations,$annotation_list);
     if($format=="json")echo json_encode(array('organisms'=>$annotation_list));
     else{echo'<organisms>'.displayXml($annotations).'</organisms>';}
 }
 elseif($genome_listing){ //the query is to list all organism with chromosome assembly data
      if(empty($host))displayOrganismWgenome($annotation_list);
      else ucsc_displayOrganismWgenome($annotation_list);$annotations=array_merge($annotations,$annotation_list);
      if($format=="json")echo json_encode(array('organisms'=>$annotation_list));
      else{echo'<organisms>'.displayXml($annotations).'</organisms>';}
  }
  elseif($annotation_load){//the query is to download a specific gene annotations file for the specified organism version
       if(empty($host))downloadAnnotations( $organism_version,$gene_prediction,$annotation_list);
       else ucsc_downloadAnnotations($organism_version,$gene_prediction,$annotation_list,$atype);
       $annotations=array_merge($annotations,$annotation_list);
       if($format=="json")echo json_encode(array('annotfiles'=>$annotation_list));
       else{echo'<annotfiles>'.displayXml($annotations).'</annotfiles>';}
  }
  elseif($chromosome_listing){
      getChromosomeList($org_vid,$annotation_list);
      if($format=="json")echo json_encode(array('chromosomes'=>$annotation_list));
      else{echo'<chromosomes>'.displayXml($annotation_list).'</chromosomes>';}
  }
  elseif(!empty($term)){//query for annotations using a term (gene name,gene id, transcript,...)
      $hit_list=array(); 
      list($chromosome_id,$chrom_start,$chrom_end)=explode(":",getGenomicRegion($term));
      if($chromosome_id>0){$hit_list=getCoordHits($org_vid,$chromosome_id,$chrom_start,$chrom_end,$strand);}
      else{
         if(preg_match("/mgi:/i",$term))$hit_list=getMgiHits($term,$org_vid);
         else $hit_list=getTermHits($term,$org_vid);
      }
      $pred_id= ($gene_prediction !="All")?implode(",",getPredictionIdList($gene_prediction)):0; 
      $annotations=array();$annotat=array(); $pg="Gene Prediction Source= $gene_prediction";$header="$head$pg";
      $annotations[]=array("head"=>$header);
      if(count($hit_list)>0){
        if($isoform)$annotat=getIsoformAnnotation($org_vid,$hit_list,$pred_id);
        else $annotat=getGenomeAnnotation($org_vid,$hit_list,$pred_id); 
        $annotations=array_merge($annotations,$annotat);
        if($format=="json"){echo json_encode(array('annotations'=>$annotat));}
        else{echo"<annotations>".displayXml($annotations)."</annotations>";}  
      }else{
       if($format=="json"){echo json_encode(array('annotations'=>$annotat));}
       else{echo"<annotations>".displayXml($annotations)."</annotations>";} 
     }
  }
  else{//display help 
   if($format=="json")echo json_encode($usage);else echo $usage;
  }
  
 ?>

