<?php

/****************************************************************************
  Author: Lucie Ndzana Hutchins
  Date : October 2013
  
  Note: This utility for php function calls


*****************************************************************************/


$dbname="";
$dir=getcwd();
$version_file="$dir/ucsc_organism_version.txt";$ugenome_list_file="$dir/ucsc_organisms_w_chromosome_data.txt";
$annot_file="$dir/ucsc_geneannotations.txt"; $loadscript="$dir/ucsc_generateConfigs.pl";
$features_file="$dir/ucsc_features.txt";
$annot_loadscript="$dir/toolsGetCordinates.pl";$perl="/usr/bin/perl";
$service_url="http://demon.jax.org/transcriptdb/web-services/webservices_uploads";
function checkConfig(){ //I will set this to be a cron job
 global $ugenome_list_file; global $loadscript;global $dir; 
 if(!file_exists("$ugenome_list_file")){exec("$perl $loadscript $dir");}
 $Diff = (time() - filectime("$ugenome_list_file"))/60/60/24;
 if ($Diff > 21) exec("$perl $loadscript $dir");//also check the last mod date < 21 days
}
function ucsc_snp($organism_version){
   $dir="$service_url/";
   if(!is_dir($dir));
}
/********* downloads the specified gene annotations from ucsc for a given organism version
 and stores the file on the server under /tmp/webservices_uploads/
***********************************************************************************/
function ucsc_downloadAnnotations($organism_version,$gene_prediction,&$annotation_list,$type){
  global $dir;global $annot_loadscript; $file="ucsc_$organism_version-$gene_prediction-transcripts.txt"; 
  if($type=="feature")$file="ucsc_$organism_version-$gene_prediction.txt";
  $annot_file="/tmp/webservices_uploads/$file";$zipfile="$annot_file.gz";global $service_url;
  if(!file_exists("$zipfile")){ 
     exec("$perl $annot_loadscript -f $annot_file -v $organism_version -a $gene_prediction -d ucsc",$output);
   }
  $Diff = (time() - filectime("$zipfile"))/60/60/24;//also check the last mod date < 21 days
  if ($Diff > 21) exec("$perl $annot_loadscript -f $annot_file -v $organism_version -a $gene_prediction -d ucsc",$output);
  $annotation_list[]=array('annotfile'=>array('name'=>"$service_url/$file.gz"));
}
/*********************************************************
 This function returns a list of organism versions 
 that have chromosome data - parses ucsc_organisms_w_chromosome_data.txt file
 check if ucsc_organisms_w_chromosome_data.txt does not exists
 it calls loadConfigs()
**********************************************************/
function ucsc_getOrgWithGenome(){
  $versions="";global $ugenome_list_file; checkConfig();
  $lines=file($ugenome_list_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if($versions=="")$versions=$version;else $versions.=",$version";
   }
 return $versions;
}
/*********************************************************
 This function returns a list of organism versions 
 that have gene annotations
**********************************************************/
function ucsc_getOrgWithAnnot(){
  $versions="";global $annot_file; checkConfig();
  $lines=file($annot_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if($versions=="")$versions=$version;else $versions.=",$version";
   }
 return $versions;
}
/************** getListGenePrediction *******************************
 get the list of annotation sources for this organism version
*********************************************************************/
function ucsc_getListGenePrediction($orgv,&$predictions){
 checkConfig(); global $annot_file;
  $lines=file($annot_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if(strtoupper($orgv)==strtoupper($version)){
         $annots=split(":",$data);
         foreach($annots as $annot)
            $predictions[]=array('prediction'=>array('name'=>$annot,'id'=>$annot));
       break;
    }
  }
}
/************** getListFeatures ********************************
 get the list of features annotation for this organism version
****************************************************************/
function ucsc_getListFeatures($orgv,&$predictions){
  checkConfig(); global $features_file; $lines=file($features_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if(strtoupper($orgv)==strtoupper($version)){
         $annots=split(":",$data);
         foreach($annots as $annot)
            $predictions[]=array('prediction'=>array('name'=>$annot,'id'=>$annot));
       break;
    }
  }
}
/********************************************************
 This function the list of organisms with gene annotations.
 The result include the organism name, and the organism version
*********************************************************/
function ucsc_displayOrganismList($db, &$annot_list){
  checkConfig(); global $annot_file; $lines=file($annot_file);
  foreach($lines as $line){
       list($orgg,$name,$version,$data)=explode(",",$line);
       if(!empty($db)){
          $annot_list[]=array('organism'=>array('name'=>$name,'version'=>$version,'organism_id'=>$name,'version_id'=>$version));
          break;
       }else $annot_list[]=array('organism'=>array('name'=>$name,'version'=>$version,'organism_id'=>$name,'version_id'=>$version));
   }  
}
/********************************************************
 This function displays the list of organism versions 
 that have genome data (chromosome sequence)
*********************************************************/
function ucsc_displayOrganismWgenome(&$annot_list){
  global $ugenome_list_file; $lines=file($ugenome_list_file);
  foreach($lines as $line){
     list($orgg,$org,$version,$data)=explode(",",$line);
     if($version!="version")$annot_list[]=array('organism'=>array('name'=>$org,'version'=>$version));
   }
}

?>
