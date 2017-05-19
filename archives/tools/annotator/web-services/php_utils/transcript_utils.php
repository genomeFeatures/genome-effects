<?php

/****************************************************************************
  Author: Lucie Ndzana Hutchins
  Date : July 2012
  
  Note: This utility for php function calls

*****************************************************************************/
require_once("user.php");require_once("global.php");

$dbname="graber_transcriptdb";

/**************************************************************
This function filter the user input for security
***************************************************************/
function filterQuery($query){
  $query=htmlspecialchars(stripslashes($query));
  $query=str_ireplace("script", "blocked", $query);
  $query=mysql_escape_string($query);

 return $query;
}
function getGenomicRegion($term){
  list($chrom,$posData)=split(":",$term);$chromosome_id=0;$chrom_start=0;$chrom_end=0;
  if($posData){$posData=str_replace(" ","",$posData);
     $chrom=str_replace("chr","",$chrom);$chrom=str_replace("CHR","",$chrom);preg_match("/(\d+.*)-(\d+.*)/",$posData,$matches);
    if(count($matches)>2){
       $chrom_start=$matches[1];$chrom_end=$matches[2];$chromosome_id=getChromosomeId($chrom);
       if($chromosome_id>0){ $mb=1000000;
          if(($chrom_start>0)&&($chrom_end>0)){ #convert cordinates into BP if needed
            if(preg_match("/\d+\.\d+/",$chrom_end))$chrom_end *=$mb ;
            if(preg_match("/\d+\.\d+/",$chrom_start))$chrom_start*=$mb;
          }
       }
     }
  }
  return "$chromosome_id:$chrom_start:$chrom_end";
}
/*********************************************************
 This function returns the chromosome id for
 a given chromsome name

**********************************************************/
function getChromosomeId($chrom){
 $query="select chromosome_id from chromosome where chromosome_name='$chrom'";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$chr_id=0;
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0){
       $row = mysql_fetch_assoc($result);
       $chr_id=$row["chromosome_id"];
    }
 }
 return $chr_id;
}
/*********************************************************
 This function returns a list of organism versions 
 that have chromosome data

**********************************************************/
function getOrgWithGenome(){
  $versions="";global $genome_list_file;$lines=file($genome_list_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if($versions=="")$versions=$version;
     else $versions.=",$version";
   }
 return $versions;
}
/*********************************************************
 This function returns the chromosome name for
 a given chromsome id

**********************************************************/
function getChromosomeName($chrom_id){
 $query="select chromosome_name from chromosome where chromosome_id=$chrom_id";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$chr_name="";
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0){
       $row = mysql_fetch_assoc($result);
       $chr_name=$row["chromosome_name"];
     }
  }
 return $chr_name;
}
/*********************************************************
 This function returns the local transcript id associated
 to the specified term, gene prediction id and the organism 
 version id
**********************************************************/
function getTranscriptId($term,$pred_id,$org_vid){
 $query="select transcript_id from transcript_by_annotation where transcript_name='$term' ";
 $query.=" and organism_version_id=$org_vid and gene_prediction_id in($pred_id)";
  global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$tx_id=0;
 //if(!$result)queryError($query,mysql_error());
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0){
       $row = mysql_fetch_assoc($result);
       $tx_id=$row["transcript_id"];
    }
  }
 return $tx_id;
}
/********************** getPredictionId ********
 This function returns the gene prediction id
 associated to the specified gene prediction name

*************************************************/
function getPredictionIdList($prediction_name){
 $pred_id=array();$list=explode(",",$prediction_name);
 $query="select gene_prediction_id from gene_prediction where gene_prediction_name ";
 if(!is_array($list))$query.="='$prediction_name'";
 else{ $query.=" in(";$i=0;
    while($i<count($list)){
       $annotation=$list[$i];if($i==0){$query.="'$annotation'";}else $query.=",'$annotation'";++$i;
     }$query.=")";
 }
 global $con;if(!$con)$con=getConnection();$result = mysql_query($query,$con);
 if($result){
    while( $row = mysql_fetch_assoc($result))$pred_id[]=$row["gene_prediction_id"];
  }
 return  $pred_id;
}
/********************** getPredictionId ********
 This function returns the gene prediction id
 associated to the specified gene prediction name

*************************************************/
function getPredictionId($prediction_name){
 $pred_id=0;
 $query="select gene_prediction_id from gene_prediction where gene_prediction_name ='$prediction_name'";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0){
       $row = mysql_fetch_assoc($result);
       $pred_id=$row["gene_prediction_id"];
    }
  }
 return  $pred_id;
}
/********************** getPredictionName ********
 This function returns the gene prediction name
 associated to the specified gene prediction id

*************************************************/
function getPredictionName($prediction_id){
 $pred_name="";
 $query="select gene_prediction_name from gene_prediction where gene_prediction_id =$prediction_id";
 global $con;if(!$con)$con=getConnection();$result = mysql_query($query,$con);
 if($result){$numRows =mysql_numrows($result);
   if($numRows>0){ $row = mysql_fetch_assoc($result);$pred_name=$row["gene_prediction_name"];}
 }
 return  $pred_name;
}
/************** getDefaultGenePrediction ***********
 get the list of annotation sources for this organism version
 and set the default annotation for this organism
****************************************************/
function getDefaultGenePrediction($orgv_id,$term,$table){
 $query="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p ";
 $query.="where op.organism_version_id =$orgv_id and op.gene_prediction_id=p.gene_prediction_id ";
 $query.=" and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
 if($table!=""){
    if($table=="gene_by_annotation"){
      $query="select distinct gene_prediction_name from gene_by_annotation g,gene_prediction p";
      $query .=" where g.gene_name ='$term' and g.organism_version_id=$orgv_id and g.gene_prediction_id=p.gene_prediction_id";
    }
   elseif($table=="transcript_by_annotation"){
      $query="select distinct gene_prediction_name from transcript_by_annotation g,gene_prediction p";
      $query .=" where g.transcript_name ='$term' and g.organism_version_id=$orgv_id and g.gene_prediction_id=p.gene_prediction_id";
    }
    elseif($table=="cc_founders_transcripts"){
       $query="select distinct gene_prediction_name from cc_founders_transcripts op, gene_prediction p";
       $query.="where op.transcript_name='$term' and op.gene_prediction_id=p.gene_prediction_id ";
    }
    elseif($table=="mgi_genes"){
      $query="select distinct gene_prediction_name from gene_by_annotation g, mgi_genes m,gene_prediction p";
      $query .=" where m.mgi_geneid ='$term' and m.mgi_symbol=g.gene_name ";
      $query .=" and g.organism_version_id=$orgv_id and g.gene_prediction_id=p.gene_prediction_id";
    }
   
 }
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$gene_prediction="";
 if($result){ $numRows =mysql_numrows($result);$predictions=array();
   if($numRows>0){
       while($row = mysql_fetch_assoc($result)){
           $gene_prediction=$row["gene_prediction_name"];$predictions{$gene_prediction}=1;
       }
       $ens="ensGene";$mgi="mgiGene";$b6="C57BLGene";
       $db=getOrganismVersionName($orgv_id);if($db=="mm10")$db="GRCm38-";
       if($predictions{"$db$ens"}==1){$gene_prediction="$db$ens";}
       elseif($predictions{"$db$mgi"}==1){$gene_prediction="$db$mgi";}
       elseif($predictions{"$db$b6"}==1){$gene_prediction="$db$b6";}
       elseif($predictions{"knownGene"}==1){$gene_prediction="knownGene";}
       elseif($predictions{"refGene"}==1){$gene_prediction="refGene";}
       elseif($predictions{"xenoRefGene"}==1){$gene_prediction="xenoRefGene";}
    }
  }
 return $gene_prediction;
}
/************** getTermGenePrediction ***********
 get the list of annotation sources for this term
 and for this organism version
****************************************************/
function getTermGenePrediction($term,$table_name,$orgv_id,$field_name){
 $query="select distinct gene_prediction_id from $table_name ";
 $query.="where $field_name='$term'  and organism_version_id =$orgv_id ";
 $query.=" and gene_prediction_id not in (1,2,4)";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$gene_prediction="";
 if($result){
    $numRows =mysql_numrows($result);$predictions=array();
    if($numRows>0){
        while($row = mysql_fetch_assoc($result)){
           $gene_prediction=$row["gene_prediction_id"];
           $predictions{$gene_prediction}=1;
        }
       if($predictions{$gene_prediction}==13)$gene_prediction=13;
       elseif($predictions{$gene_prediction}==3)$gene_prediction=3;
       elseif($predictions{$gene_prediction}==6)$gene_prediction=6;
     }
  }
 return $gene_prediction;
}

/************** checkTermGenePrediction ***********
 check the current term was provided by this annotation source
 and for this organism version
****************************************************/
function checkTermGenePrediction($term,$table_name,$orgv_id,$field_name,$pred_id){
 $query="select * from $table_name ";
 $query.="where $field_name='$term' and organism_version_id =$orgv_id ";
 $query.=" and gene_prediction_id in($pred_id)";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$gene_prediction=0;
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0)$gene_prediction=1;
  }
 return $gene_prediction;
}
/************** getListGenePrediction ***********
 get the list of annotation sources for this organism version
****************************************************/
function getListGenePrediction($orgv_id,&$predictions){
 $query="select distinct gene_prediction_name as name,p.gene_prediction_id";
 $query.=" from gene_prediction_by_organism op, gene_prediction p "; //gene_prediction_by_organism
 $query.=" where op.organism_version_id =$orgv_id and op.gene_prediction_id=p.gene_prediction_id ";
 $query.=" and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
 global $con;if(!$con)$con=getConnection();$result = mysql_query($query,$con);
 //if(!$result)queryError($query,mysql_error());
 if($result){
    while($row = mysql_fetch_assoc($result)){
      $predictions[]=array('prediction'=>array('name'=>$row["name"],'id'=>$row["gene_prediction_id"]));
    }
  }
}
/****************** getChromosomeList *************
**************************************************/
function getChromosomeList($org_vid,&$chromosomes){
 $query="select chromosome_name as name, c.chromosome_id as id,chrom_size as size ";
 $query.=" from chromosome c, chromosome_data d where d.organism_version_id=$org_vid ";
 $query.=" and d.chromosome_id=c.chromosome_id";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);
 if($result){
    while($row = mysql_fetch_assoc($result)){
          $chromosomes[]=array('chromosome'=>array('name'=>$row["name"],'id'=>$row["id"],'size'=>$row["size"]));
    }
 }
}
/*************** getOrganism ***********************
 Given the local organism version id, this function returns 
 a ucsc database name
*******************************************************/
function getOrganismVersionName($org_vid){
  global $con;if(!$con)$con=getConnection();
  $query="select ucsc_db from organism_version where organism_version_id=?";
  $result = mysql_query($query,$con);$org_version_id=0;
  if($result){
     $numRows =mysql_numrows($result);
     if($numRows>0){
        $row = mysql_fetch_assoc($result);
        $ucsc_db=$row["ucsc_db"];
     }
  }
 return $ucsc_db;
}
/*************** getOrganismVid ***********************
 Given a ucsc database name, this function returns 
the local organism version id
*******************************************************/
function getOrganismVid($ucsc_db){
  global $con;global $dbname;if(!$con)$con=getConnection();
  $query="select organism_version_id from organism_version where ucsc_db='$ucsc_db'";
  $result = mysql_query($query,$con);$org_version_id=0;
  if($result){
     $numRows =mysql_numrows($result);
     if($numRows>0){
        $row = mysql_fetch_assoc($result);
        $org_version_id=$row["organism_version_id"];
     }
  }
 return $org_version_id;
}
/**** getItemsType    *********************************************
***** This function returns the type of a given search term *******
***** graber_transcript has the following types:            *******
      chromosome_region -> 1
      gene id ->  2
      transcript id ->  3
      EST -> 4
******************************************************************/
function getItemsType($searchTerm){
   $searchTerm =trim($searchTerm);global $con;
   if(!$con)$con=getConnection();
   $query="select distinct term as term_name,term_type_id as term_type from query_terms";
   $query .=" where  term like '$searchTerm%' ";
   $result = mysql_query($query,$con);$itemTypeList=array();
   if($result){ 
      $numRows =mysql_numrows($result);
      while($row = mysql_fetch_assoc($result)){
            $term=$row["term_name"];$id=$row["term_type"];$itemTypeList[$id]=$term;
      }
    }
   return $itemTypeList;
}
function getMgiHits($searchTerm,$organism_version_id){
$query=" select distinct t.transcript_id,transcript_name,t.gene_prediction_id ";
$query.=" from gene_by_annotation g, mgi_genes m, transcript_by_annotation t ";
$query .=" where m.mgi_id ='$searchTerm' and m.mgi_symbol=g.gene_name and g.organism_version_id=$organism_version_id";
$query.=" and g.transcript_id=t.transcript_id and  g.organism_version_id=t.organism_version_id ";

$result = mysql_query($query,$con);$hits=array();
if($result){
   while($row = mysql_fetch_assoc($result)){
      $hits[]=array('id'=>$row["transcript_id"],
                    'name'=>$row["transcript_name"],
                    'prediction_id'=>$row["gene_prediction_id"]);
    }
 }
 return $hits;
}
/**** getTermHits    *********************************************
***** This function returns the hits of a given search term *******
***** graber_transcript has the following types:            *******
   
******************************************************************/
function getTermHits($searchTerm,$organism_version_id){
   $searchTerm =trim($searchTerm);global $con;$hits=array();
   if(!$con)$con=getConnection();
   $query="select distinct t.transcript_id,transcript_name,t.gene_prediction_id   from gene_by_annotation g, transcript_by_annotation t 
        where gene_name ='$searchTerm' and g.organism_version_id=$organism_version_id and g.transcript_id=t.transcript_id
        and  g.organism_version_id=t.organism_version_id";
   $result = mysql_query($query,$con);
   if($result){
      while($row = mysql_fetch_assoc($result)){
           $hits[]=array('id'=>$row["transcript_id"],
                        'name'=>$row["transcript_name"],
                         'prediction_id'=>$row["gene_prediction_id"]);
      }
    }
   if(count($hits)<=0){ 
     $query="select distinct transcript_id,transcript_name,gene_prediction_id from transcript_by_annotation ";
      $query .=" where transcript_name ='$searchTerm' and organism_version_id=$organism_version_id  ";
      $result = mysql_query($query,$con);
      if($result){
         while($row = mysql_fetch_assoc($result)){
                $hits[]=array('id'=>$row["transcript_id"],
                        'name'=>$row["transcript_name"],'prediction_id'=>$row["gene_prediction_id"]);
           }
       } 
   }
   elseif(count($hits)<=0){  //check if cc founders transcript_name
      $query="select distinct t.transcript_id,t.transcript_name,t.gene_prediction_id 
        from transcript_by_annotation t,cc_founders_transcripts c
        where  c.transcript_name='$searchTerm' and t.organism_version_id=$organism_version_id  and c.gene_prediction_id = t.gene_prediction_id 
        and    t.transcript_name=c.transcript_id ";
      $result = mysql_query($query,$con);
      if($result){ 
         while($row = mysql_fetch_assoc($result)){ 
                $hits[]=array('id'=>$row["transcript_id"],
                        'name'=>$row["transcript_name"],
                         'prediction_id'=>$row["gene_prediction_id"]);
         }
      }
   }
   elseif(count($hits)<=0){  //check if cc founders gene_id
      $query="select distinct t.transcript_id,t.transcript_name,t.gene_prediction_id"; 
      $query.=" from gene_by_annotation g,cc_founders_genes c,transcript_by_annotation t
        where  c.gene_id ='$searchTerm' and c.gene_name=g.gene_name and g.gene_prediction_id=c.gene_prediction_id 
        and g.organism_version_id=$organism_version_id  and g.transcript_id=t.transcript_id 
        and  g.organism_version_id=t.organism_version_id "; 
      $result = mysql_query($query,$con);
      if($result){ 
         while($row = mysql_fetch_assoc($result)){ 
                $hits[]=array('id'=>$row["transcript_id"],
                        'name'=>$row["transcript_name"],
                         'prediction_id'=>$row["gene_prediction_id"]);
         }
      }
   }
   return $hits;
}
/***************************************************************************
 This function returns the transcript ids of all the transctipts 
 that overlap the specified genomic range  c.chromosome_name,t.transcript_id,tx_start,tx_end,strand 
****************************************************************************/
function getCoordHits($organism_version_id,$chromosome_id,$tStart,$tEnd,$strand){
  global $con;
  if(!$con)$con=getConnection();$query="";$hits=array();
  $query=" select  distinct c.transcript_id,t.transcript_name,t.gene_prediction_id ";
  $query.=" from transcript c,transcript_by_annotation t";
  $query.=" where c.organism_version_id=$organism_version_id ";
  $query.=" and c.chromosome_id=$chromosome_id  ";
  if($strand!="")$query.=" and strand='$strand'";
  $query.=" and tx_end >= $tStart and tx_start<=$tEnd ";
  $query.=" and t.transcript_id=c.transcript_id limit 500";
  $result = mysql_query($query,$con);
  if(!$result)queryError($query,mysql_error());
  if($result){
     while($row = mysql_fetch_assoc($result)){
           $hits[]=array('id'=>$row["transcript_id"],'name'=>$row["transcript_name"],
                         'prediction_id'=>$row["gene_prediction_id"]);
     } mysql_free_result($result);
   }
  
  return $hits;
}

/***************************************************************************
 This function returns the cordinates of all the transctipts 
 for the specified transcript id
****************************************************************************/
function getOverLapTranscripts($organism_version_id,$tx_list){ 
  global $con;if(!$con)$con=getConnection();$query="";
  $query=" select c.chromosome_name,c.chromosome_id,t.transcript_id,tx_start,tx_end,strand ";
  $query.=" from transcript t,chromosome c where organism_version_id=$organism_version_id ";
  $query.=" and transcript_id in($tx_list) and t.chromosome_id=c.chromosome_id ";
  $query.="order by t.chromosome_id,tx_start,t.transcript_id";
  $result = mysql_query($query,$con);$data=array(); if(!$result)queryError($query,mysql_error());
  if($result){
     while($row = mysql_fetch_assoc($result)){
           $data[]=array('tid'=>$row["transcript_id"],
                         'tStart'=>$row["tx_start"],
                         'tEnd'=>$row["tx_end"], 'chrom'=>$row["chromosome_name"],
                         'chrid'=>$row["chromosome_id"],
                         'strand'=>$row["strand"]);
     } mysql_free_result($result);
   }
  
  return $data;
}
/***************************************************************************
 This function returns MGI gene symbol of a given MGI ID 
****************************************************************************/
function getMGIgene($mgiID){
  $query="select mgi_symbol from mgi_genes where mgi_id='$mgiID'";
  global $con;if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);$data="";
  if($result){
     while($row = mysql_fetch_assoc($result)){$data=$row["mgi_symbol"];}
     mysql_free_result($result);
  }
  return $data;
}
/***************************************************************************
 This function returns gene symbols provided by a give annotation source
 for a given transcriptID 
****************************************************************************/
function getTxGene($tx_id,$pred_id){
  $query="select distinct gene_name from gene_by_annotation "; $query.=" where transcript_id=$tx_id ";
  if($pred_id>0)$query.=" and gene_prediction_id in($pred_id)";
  global $con;if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);$data=array();
  if($result){
     while($row = mysql_fetch_assoc($result)){
            $data[]=array('gene'=>$row["gene_name"]);
     }mysql_free_result($result);return $data;
  }
}
/***************************************************************************
 This function returns all the exon coordinates of a given transcriptID and 
  gene prediction source is specified 
****************************************************************************/
function getTxExons($tx_id,$tx,$pred_id,$org_vid){
   $query="select distinct e.exon_id,chromosome_name,exon_start,exon_end,strand,exon_frame from transcript_exon t,exon e,chromosome c ";
   $query.=" where  transcript_id in($tx_id) and t.organism_version_id=$org_vid and gene_prediction_id in ($pred_id)";
   if(!empty($tx))$query.=" and transcript_name='$tx'";
   $query.=" and t.organism_version_id=e.organism_version_id and t.exon_id=e.exon_id and e.chromosome_id=c.chromosome_id";
   $query.=" order by chromosome_name,strand,exon_start asc";
   global $con;if(!$con)$con=getConnection();
    $result = mysql_query($query,$con);$data=array();
   if($result){
     while($row = mysql_fetch_assoc($result)){
           $data[]=array("id"=>$row["exon_id"],
                        'chrom'=>$row["chromosome_name"],'ex_start'=>$row["exon_start"],
                         'ex_end'=>$row["exon_end"],'ex_frame'=>$row["exon_frame"],'strand'=>$row["strand"]);
      }mysql_free_result($result);
   }
 return $data;
}
/***************************************************************************
 This function returns all the exons that overlap a given genomic region and 
 if gene prediction source is specified 
****************************************************************************/
function getOverlapExons($strand,$chromosome_id,$tStart,$tEnd,$organism_version_id,$transcriptList){
  $query=" select distinct chromosome_id,e.exon_id,e.exon_start,e.exon_end,e.strand ";
  $query.=" from exon e, transcript_exon te where e.organism_version_id=$organism_version_id ";
  $query.=" and chromosome_id=$chromosome_id and strand='$strand'";
  $query.=" and exon_end >=$tStart and exon_start<=$tEnd ";
  $query.=" and e.organism_version_id=te.organism_version_id and e.exon_id=te.exon_id";
  $query.=" and te.transcript_id in($transcriptList) order by exon_start";
  global $con;if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);$data=array();
  //if(!$result)queryError($query,mysql_error());
  if($result){
     while($row = mysql_fetch_assoc($result)){
            $data[]=array("id"=>$row["exon_id"],
                   'chrom'=> $row["chromosome_id"],
                   'tStart'=>$row["exon_start"],
                   'tEnd'=>$row["exon_end"],
                   'strand'=>$row["strand"]);
     }mysql_free_result($result);
  }
  return $data;
}

/****************************************************************
 This function returns all the cds starts and ends of a given 
 transcript from the provided gene predictions list
*****************************************************************/
function getCDS($tx_name,$tx_id,$orgvid,$predict_id){
 $query="select distinct cdsStart,cdsEnd from transcript_translation ";
 $query.=" where organism_version_id=$orgvid and transcript_name='$tx_name'";
 $query.= " and transcript_id=$tx_id and gene_prediction_id in($predict_id)";
 global $con;if(!$con)$con=getConnection();
 $cds=array();$cdsStarts;$cdsEnds;
  $result = mysql_query($query,$con);$data=array();
  if(!$result)queryError($query,mysql_error());
  if($result){
      while($row = mysql_fetch_assoc($result)){
        $cds_start=$row["cdsStart"]+1;
        if(empty($cdsStarts)){$cdsStarts=$cds_start;$cdsEnds=$row["cdsEnd"];}
        else {$cdsStarts.=",".$cds_start;$cdsEnds.=",".$row["cdsEnd"];}
      }
     $cds["start"]=$cdsStarts;$cds["end"]=$cdsEnds;
  }
 return $cds;
}
/***************************************************************************
 This function returns the isoform view of the query result
****************************************************************************/
function getIsoformAnnotation($org_vid,$tx_list,$pred_id){
   $annotation=array();$region_start=0;$region_end=0; $txmap=array();$genes=array();$tx="";
   while(list($index,$row)=each($tx_list)){
       $tx_id=$row["id"];$name=$row["name"];$prediction_id=$row["prediction_id"];
       if(($prediction_id<3))continue;
       if(($pred_id<=0)&&($prediction_id==5))continue;//skip xenoref  by default
       if(($pred_id<=0)&&($prediction_id>24&&$prediction_id<=40))continue;//skip cc founders by default
       $row["tx_start"]+=1; if($pred_id>0){if($pred_id!=$prediction_id)continue;}
       $pred_name=getPredictionName($prediction_id);$tx_coord=getOverLapTranscripts($org_vid,$tx_id);
       $tx_coord[0]["tStart"]+=1;$gene_names=""; $ex_starts=""; $ex_ends="";$ex_frames="";
       $strand= $tx_coord[0]["strand"];$cds=getCDS($name,$tx_id,$org_vid,$prediction_id);
       $qh_getTxGene=getTxGene($tx_id,$prediction_id);
       while(list($index2,$gene)=each($qh_getTxGene)){ $gene_symbol=$gene["gene"];
          if($gene_names == ""){$gene_names.="$gene_symbol";}
          else{$gene_names.=",$gene_symbol";}
       }$ex_start=0;$ex_end=0;$qh_getTxExon=getTxExons($tx_id,$name,$prediction_id,$org_vid);
       while(list($index3,$value)=each($qh_getTxExon)){
             $ex_start=$value["ex_start"]+1;$ex_end=$value["ex_end"];$ex_frame=$value["ex_frame"];
             if($ex_starts ==""){$ex_starts="$ex_start";$ex_ends="$ex_end";$ex_frames="$ex_frame";}
             else{$ex_starts.=",$ex_start";$ex_ends.=",$ex_end";$ex_frames.=",$ex_frame";}
        }
        if($pred_id>24&&$pred_id<=40)$ex_frames="";
        $annotation[]=array("annotation"=>array("source"=>$pred_name,"chrom"=>$tx_coord[0]["chrom"],
                          "txStart"=>$tx_coord[0]["tStart"],"txEnd"=>$tx_coord[0]["tEnd"],
                          "strand"=>$strand,"Name"=>$name, "Name2"=>$gene_names,
                          "cdsStart"=>$cds["start"],"cdsEnd"=>$cds["end"],
                          "exonStarts"=>$ex_starts,"exonEnds"=>$ex_ends,"exonFrames"=>$ex_frames));
  }
  return $annotation;
}
/***************************************************************************
 This function returns annotations that overlap the specified genomic range
 in aggregate format
****************************************************************************/
function getGenomeAnnotation($org_vid,$tx_list,$pred_id){ 
  $region_map=array();$tx_map=array();$getTxAccession=array();$annotation=array();$qh_getTxGene=array();
  $geneAccession=array();$pred_list=array();
  while(list($index,$row)=each($tx_list)){ //generate a list of tx ids
       if(($row["prediction_id"]<3))continue;    //skip all_est, and all_mrna
       if(($pred_id<=0)&&($row["prediction_id"]==5))continue;   //skip xenoref  by default
       if(($pred_id<=0)&&($row["prediction_id"]>24&&$row["prediction_id"]<=40))continue;//skip cc founders by default
       if(($pred_id!=0)&&($pred_id!=$row["prediction_id"]))continue; //only display selected annotations
       $pred_list{getPredictionName($row["prediction_id"])}=1;
       $getTxAccession{$row["id"]}{$row["name"]}=1;$qh_getTxGene=getTxGene($row["id"],$row["prediction_id"]);
       while(list($index2,$gene)=each($qh_getTxGene)){$geneAccession{$row["id"]}{$gene["gene"]}=1;}
  }$qh_Overlaplist=getOverLapTranscripts($org_vid,implode(",",array_keys($getTxAccession)));
  $current_start=0;$current_end=0; $ntx_list=""; $chromosome="";$strand="";$chrid=0;
  $pred_name=implode(",",array_keys($pred_list));
  //generate unique genomic regions -- aggregate transcripts
  while(list($index,$value)=each($qh_Overlaplist)){ 
     $tx_start=$value["tStart"];++$tx_start;$tx_end=$value["tEnd"];$tx_id=$value["tid"];
     $chromosome=$value["chrom"]; $strand=$value["strand"];$chrid=$value["chrid"];
     if($current_start==0){$current_start=$tx_start;$current_end=$tx_end;$tx_map=array();$tx_map{$tx_id}=1;}
     else{  //transcripts are stored on the map sorted by txStart asc
         if(($tx_start>=$current_start)&&($tx_start<=$current_end)){
            if($current_end<$tx_end){$current_end=$tx_end;}
            if(!array_key_exists($tx_id,$tx_map)){$tx_map{$tx_id}=1;}
         }else{ #this is a new genomic region, there is no overlap with previous transcript region
            $region_map{"$chromosome:$chrid"}{$strand}{$current_start}{"cord"}=$current_end;
            $region_map{"$chromosome:$chrid"}{$strand}{$current_start}{"tx"}=implode(",",array_keys($tx_map));
            $current_start=$tx_start;$current_end=$tx_end;$tx_map=array();$tx_map{$tx_id}=1;
         }
     }
  }//store the last region
  $region_map{"$chromosome:$chrid"}{$strand}{$current_start}{"cord"}=$current_end;
  $region_map{"$chromosome:$chrid"}{$strand}{$current_start}{"tx"}=implode(",",array_keys($tx_map));
  #now loop through the unique transcript list to get the unique exons
  ksort($region_map); 
  foreach($region_map as $chromosome=>$strandval){ksort($strandval);
          foreach($strandval as $strand => $region_startVal){ ksort($region_startVal);
             foreach($region_startVal as $txstart=>$row){$transcriptList="";$geneList="";
                     $txend=$row{"cord"}; $transcriptLists=array();$geneLists=array();$tx_lists=explode(",",$row{"tx"});
                  //generate transcript accessions and genes
                  foreach($tx_lists as $tx_id){
                     $keys=array_keys($getTxAccession{$tx_id});foreach($keys as $tx){$transcriptLists{strtoupper($tx)}=1;}
                     $keys=array_keys($geneAccession{$tx_id});foreach($keys as $gene){$geneLists{strtoupper($gene)}=1;}
                  } //get aggregate exons for this region
                  $transcriptList=implode(",",array_unique(array_keys($transcriptLists)));
                  $geneList=implode(",",array_unique(array_keys($geneLists)));
                  $ex_starts="--"; $ex_ends="--";list($chrom,$chrom_id)=explode(":",$chromosome);
                  $current_start=0;$current_end=0;$ex_start=0;$ex_end=0;
                  $qh_getOverlapE=getOverlapExons($strand,$chrom_id,$txstart,$txend,$org_vid,$row{"tx"});
                  while(list($index,$value)=each($qh_getOverlapE)){
                      $ex_start=$value["tStart"];++$ex_start;$ex_end=$value["tEnd"];$exon_id=$value["id"];
                      if($current_start==0){$current_start=$ex_start;$current_end=$ex_end;}
                      else{
                          if(($ex_start>=$current_start)&&($ex_start<=$current_end)){
                             if($current_end<$ex_end){$current_end=$ex_end;}
                          }
                         else{ #this is a new genomic region, there is no overlap with previous exon
                             if($ex_starts =="--"){ $ex_starts="$current_start";$ex_ends="$current_end";}
                             else{$ex_starts.=",$current_start";$ex_ends.=",$current_end";}
                             $current_start=$ex_start;$current_end=$ex_end;
                         }
                       }
                 } #now add the last exon
                if($current_start==$ex_start){
                   if($ex_starts =="--"){$ex_starts="-$current_start-";$ex_ends="$current_end";}
                   else{$ex_starts.=",$current_start";$ex_ends.=",$current_end";}
                }
               //the aggregate version won't have cds info
               $annotation[]=array("annotation"=>array("source"=>$pred_name,"chrom"=>$chrom,
                          "txStart"=>$txstart,"txEnd"=>$txend,
                          "strand"=>$strand,"Name"=>$transcriptList, "Name2"=>$geneList,
                          "exonStarts"=>$ex_starts,"exonEnds"=>$ex_ends));
                 ////
            } // end of  foreach($region_startVal as $txstart=>$row){
         }//foreach($strandval as $strand => $region_startVal){
  }
 return $annotation;
}
/********************************************************
 This function the list of organisms with gene annotations.
 The result include the organism name, and the organism version
*********************************************************/
function displayOrganismList($db, &$annot_list){
   $query="select distinct o.organism_name as name,ucsc_db as version,";
   $query.=" o.organism_id,ov.organism_version_id as version_id";
   $query.=" from organism o,organism_version ov";
   $query.=" where o.organism_id=ov.organism_id ";
   if(!empty($db)){
       $list=array(); $list=explode(",",$db);$query.=" and ucsc_db in(";$i=0;
       while($i<count($list)){$db=$list[$i];
             if($i==0)$query.="'$db'";else $query.=",'$db'";++$i;
       }$query.=")";
    }$query.=" order by organism_name";global $con;
   if(!$con)$con=getConnection();$result =mysql_query("$query",$con);$data= array();
   if($result){
      while($row = mysql_fetch_assoc($result))
      {   $name=$row["name"];$version=$row["version"];$org_id=$row["organism_id"];$version_id=$row["version_id"];
          $annot_list[]=array('organism'=>array('name'=>$name,'version'=>$version,'organism_id'=>$org_id,'version_id'=>$version_id));
       }mysql_free_result($result);
   }  
}
/********************************************************
 This function displays the list of organism versions 
 that have genome data (chromosome sequence)
*********************************************************/
function displayOrganismWgenome(&$annot_list){
  global $genome_list_file; $lines=file($genome_list_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if($version!="version")$annot_list[]=array('organism'=>array('name'=>$org,'version'=>$version));
   }
}
function displayHelp($file){
 // set up the document
$xml = new XmlWriter();
$xml->openMemory();
//$xml->startDocument('1.0', 'UTF-8');
$xml->startElement('doc');

// CData output
$xml->startElement('help');
$xml->writeCData("$file");
$xml->endElement();

// end the document and output
$xml->endElement();
echo $xml->outputMemory(true);
}
/**********************************************************
 This function displays the xml version of the annotations
 an associative array
***********************************************************/
function displayXml($annotations){
 $data="";
  foreach($annotations as $index => $annotation) {
    if(is_array($annotation)) {
       foreach($annotation as $key => $value) {
           $data.="<$key>";
           if(is_array($value)) {
              foreach($value as $tag => $val) {
                  if(is_array($val)){
                     foreach($val as $tag => $vals){
                        $data.="<$tag>$vals</$tag>\n";
                     }
                  }
                  else{$data.="<$tag>$val</$tag>\n";}
               }
            }
            else{$data.= $value;}
            $data.="</$key>\n\n";
        }
      }
   }
 return $data;
}

?>
