<?php 
 
 $filename="snp-service.xml";
 $format = strtolower($_GET['format']) == 'json' ? 'json' : 'xml'; #xml format is default
 if($format=="json")$filename="snps-service.json";
// header("Content-Disposition: attachment; filename=\"$filename\"");
 $chrom_start=isset($_GET['chromStart'])? intval($_GET['chromStart']) : 0;
 $chrom_end=isset($_GET['chromEnd'])? intval($_GET['chromEnd']) : 0;
 $source_listing=isset($_GET['src'])? intval($_GET['src']) : 0; 
 $strain_listing=isset($_GET['str'])? intval($_GET['str']) : 0; 
 $source_id=isset($_GET['source'])? intval($_GET['source']) : 21;         // default to Imputed - Jeremy R. Wang et al. 2012
 $function_class=isset($_GET["fclass"])? $_GET["fclass"] : ""; //default functional SNPs 
 $geno=isset($_GET["geno"])? intval($_GET["geno"]) : 0; //default genotype
 $strand=$_GET['strand']; 
 header("Content-type: application/octet-stream");header("Cache-Control: cache, must-revalidate");
 header("Pragma: public");header("Content-type: application/$format");
 require_once("cgd_utils.php");
 $term=isset($_GET["t"])?filterQuery($_GET["t"]):"";$mb=1000000;$gene_id=0;$tx_id=0;
 $chrom=isset($_GET["c"])?filterQuery(strtoupper($_GET["c"])):"";$pattern = '/CHR/i';
 $snplist=array();
 if(!empty($term)){ $gene_id=getGeneID($term);$tx_id=getTranscriptID($term);
    if(($gene_id==0)&&($tx_id==0)){   //get chromosome data
       list($chrom,$posData)=explode(":",$term);list($chrom_start,$chrom_end)=explode("-",$posData);
    }
 }
 if(($gene_id>0)||($tx_id>0)){$chrom_start=0;
    getTxSNPs($tx_id,$gene_id,$function_class,$source_id,$snplist,$geno);
 }
 else if(!empty($chrom)){$chrom=preg_replace($pattern,"",$chrom);
    $chr_id=getChromosomeId($chrom);
    if($chr_id>0){
       if(($chrom_start>0)&&($chrom_end>0)){
         if(preg_match("/\d+\.\d+/",$chrom_end))$chrom_end *=$mb ;#convert cordinates into BP if entered as MBp
         if(preg_match("/\d+\.\d+/",$chrom_start))$chrom_start*=$mb;$gene_id=0;$tx_id=0;
         $txIds=getTranscriptIDsOverlap($chr_id,$chrom_start,$chrom_end,$strand);$gene_id=0;
         if(count( $txIds)>0){ 
           foreach($txIds as $tx_id){
             getTxSNPs($tx_id,$gene_id,$function_class,$source_id,$snplist,$geno);} 
        }
      }
    }
  }
 if(count($snplist)>0){
    $snps['snps']= $snplist;
    if($format=="json"){echo json_encode($snps);}
    else{
       $xml=new SimpleXMLElement("<?xml version='1.0'?><snps></snps>");
       echo"<snps>".displayXml($snplist)."</snps>"; 
     //array_to_xml($snplist,$xml);echo $xml->asXML();
   } 
 }
 elseif($source_listing){
     $sourcelist=getSourceList();
     if($format=="json")echo json_encode(array('sources'=>$sourcelist));
     else{echo'<sources>'.displayXml($sourcelist["sources"]).'</sources>';}
  }
  elseif($strain_listing){
     $anotations=getSourceStrainList($source_id);
     if($format=="json")echo json_encode(array('strains'=>$anotations));
     else{echo'<strains>'.displayXml($anotations["strains"]).'</strains>';}
  }
  else{
     $snps['snps']= $snplist;
    if($format=="json"){echo json_encode($snps);}
    else{
       $xml=new SimpleXMLElement("<?xml version='1.0'?><snps></snps>");
       echo"<snps>".displayXml($snplist)."</snps>"; 
   } 
 }
 
 ?>

