<?php

/****************************************************************************
  Author: Lucie Ndzana Hutchins
  Date : July 2008
  
  Note: This utility for php function calls

*****************************************************************************/
require_once("user.php");
require_once("config.php");

$snp_service="http://demon.jax.org/transcriptdb/web-services/snps.php";
$txdb_service_json="http://demon.jax.org/transcriptdb/web-services/index.php?format=json";
$txdb_service_xml="http://demon.jax.org/transcriptdb/web-services/index.php?format=xml";

/**************************************************************
This function filter the user input for security
***************************************************************/
function filterQuery($query){
  $query=htmlspecialchars(stripslashes($query));
  $query=str_ireplace("script", "blocked", $query);
  $query=mysql_escape_string($query);
 return $query;
}
/******************************************************
 This function displays the search section of the page
*******************************************************/
function getToolDisplay($organism_version,$annotations,$query){
 global $queryTypes; if(empty($query))$query="2:105509053-105543317";
 $data='<div id="tool-main"><div class="tool" style="z-index=0;"><table><tr>';
  $data.='<td style="padding-left:3em;padding-top:2em;padding-bottom:1em;">
            <span><h3>Enter A Genomic region of interest Or a Term
              <a name="qhint" class="qhint" onmouseover="javascript:document.getElementById(\'qhint\').style.display=\'block\';" 
                          onmouseout="document.getElementById(\'qhint\').style.display=\'none\';"> (?) </a></h3></span><br/>
            <div id="qhint"><h4>You can query using the following terms:</h4>
                           <span><b>A gene symbol or id: </b> (pax6,Cnksr3,ENSMUSG00000015202,MGI:3801900,...)</span>
                           <br/><span><b>A trancript symbol or id :</b> (ENSMUST00000144335, ...) </span>
                           <br/><span> <b>Genomic range (bp or Mbp) :</b>
                                       10: 3134304 - 3227479 OR 10: 3.1 - 3.2</span>
                            
                        </div>';
  $data.="<input type=\"txt\" name=\"q\" id=\"q\" size=\"50\" value=\"$query\"";
  $data.='onclick=\'this.value="";\'/><br/>
         </td>';
 $data.="<td>$options";
 $data.='<select name="org" id="organism" onchange="updateAnnotList();">';
 $data.='<option value="0" selected="selected">Select The Reference Organism Version</option>';
 $datas.=displayOrganismList($query_type); 
 $organisms=split(",",$datas);
  foreach($organisms as $organismdata){
   list($u_oversion,$org,$db)=split(":",$organismdata);
   if($db=="$organism_version")
       $data.="<option value=\"$db\" selected>$org &nbsp;($db )</option>";
   else $data.="<option value=\"$db\">$org &nbsp; ($db)</option>";
 }
 $data.="</select><br/><br/>";
 $data.='<input type="button" class="sb" id="qsubmit" name="submit" value="submit" onclick="validate();getFile();">
          <input type="reset" name="reset" class="sb"></div></td></tr>';
 $data.='<tr><td colspan="2"><div id="filter">';
 $data.='</div></td></tr></table></div>';
 return $data;
}
/********************************************************
 This function returns the list of current organisms 
 that have gene annotations in our database

*********************************************************/
function displayOrganismList($querytype){
   global $txdb_service_xml;$file=$txdb_service_xml."&s=1&querytype=$querytype";$data="";
   $xml = simplexml_load_file("$file");
   if($xml){ $organisms=$xml->children()->organism;$i=0; 
     while($i<count($organisms)){
          $name=$organisms[$i]->children()->name;$db=$organisms[$i]->children()->version;
           $version=$organisms[$i]->children()->version_id;
           if($i==0){$data="$version:$name:$db";}
           else $data.=",$version:$name:$db";
          ++$i;
      }
 }
 return $data; 
}
/***********************************************************************************
 This function calls the web service to extract the list of transcripts that overlap 
 a given query.The web service returns data in json format.
 The argument to the web service are:
 t=$query  => the query term
 v=$organism_version
 a=$annotationSource
************************************************************************************/
function getJsonTxListing($query,$organism_version,$annotationSource){
  global $txdb_service_json;$url="$txdb_service_json&t=$query&v=$organism_version&isoform=1";
  if(!empty($annotationSource))$url.="&a=$annotations";
  $rows="";$json_array=json_decode(file_get_contents($url),true);
  $annotationArray=$json_array{"annotations"};
  foreach($annotationArray as $index =>$data){
          $chrom=$data["annotation"]["chrom"];$tx_start=$data["annotation"]["txStart"];
          $tx_end=$data["annotation"]["txEnd"];
          $strand=$data["annotation"]["strand"]; $transcript=$data["annotation"]["Name"];
          $gene=$data["annotation"]["Name2"]; $cds_start=$data["annotation"]["cdsStart"];
          $cds_end=$data["annotation"]["cdsEnd"];$source=$data["annotation"]["source"];
          $ex_starts=$data["annotation"]["exonStarts"];$ex_end=$data["annotation"]["exonEnds"];
          $exon_count=explode(",",$ex_starts);
          $rows.="<tr $bg><td>$transcript</td><td>$chrom</td><td>$strand</td><td>$tx_start</td>";
          $rows.="<td>$tx_end</td><td>".count($exon_count)."</td>";
          $rows.="<td>$cds_start</td><td>$cds_end</td><td>$ex_starts</td><td>$ex_end</td>";
          $rows.="<td>$gene</td><td>$source</td></tr>";
  }
 return $rows; 
}
/***********************************************************************************
 This function calls the web service to extract the list of transcripts that overlap 
 a given query.The web service returns data in xml format.
 The argument to the web service are:
 t=$query  => the query term
 v=$organism_version
 a=$annotationSource
************************************************************************************/
function getXmlTxListing($query,$organism_version,$annotationSource){
  global $txdb_service_xml; 
  $url="$txdb_service_xml&t=$query&v=$organism_version&isoform=1";if(!empty($annotationSource))$url.="&a=$annotations";
  $xml = simplexml_load_file($url);$rows="";
   if($xml){
     $items= $xml->children()->annotation; $index=0; $count=count($items); $head= $xml->children()->head;
     while($index < $count){$bg="";
         if($index%2==0)$bg="style='background:#eee;'";
          $item=$items[$index]; ++$index;
          $chrom= $item->children()->chrom;$tx_start=$item->children()->txStart;$tx_end=$item->children()->txEnd;
          $strand=$item->children()->strand; $transcript=$item->children()->Name;
          $gene=$item->children()->Name2; $cds_start=$item->children()->cdsStart;
          $cds_end=$item->children()->cdsEnd;$source=$item->children()->source;
          $ex_starts=$item->children()->exonStarts;$ex_end=$item->children()->exonEnds;
          $exon_count=explode(",",$ex_starts);
          $rows.="<tr $bg><td>$transcript</td><td>$chrom</td><td>$strand</td><td>$tx_start</td>";
          $rows.="<td>$tx_end</td><td>".count($exon_count)."</td>";
          $rows.="<td>$cds_start</td><td>$cds_end</td><td>$ex_starts</td><td>$ex_end</td>";
          $rows.="<td>$gene</td><td>$source</td></tr>";
       }
    }
  return $rows; 
}
?>
