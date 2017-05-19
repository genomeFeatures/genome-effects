<?php

/****************************************************************************
  Author: Lucie Ndzana Hutchins
  Date : July 2008
  
  Note: This utility for php function calls

*****************************************************************************/
require_once("user.php");
require_once("config.php");

$dbname="graber_transcriptdb";
$xml_file="global.xml";
$snp_service="http://demon.jax.org/transcriptdb/web-services/snps.php";
$txdb_service="http://demon.jax.org/transcriptdb/web-services/index.php";
$queryTypes["a"]="Select Query Type";
//$queryTypes["txsnp"]="Overlap SNPs On Transcripts (Genes)";
$queryTypes["fsnp"]="Overlap Functional SNPs";
$queryTypes["ssnp"]="Overlap SNPs Coding Synonymous";
$queryTypes["nssnp"]="Overlap SNPs Coding NonSynonymous";
$queryTypes["tx"]="Overlap Transcripts";
//$queryTypes["txs"]="Overlap Transcript Start Sites";
$queryTypes["exon"]="Overlap Exons";
//$queryTypes["texon"]="Overlap Terminal Exons";
$queryTypes["cds"]="Overlap Coding Exons";
//$queryTypes["mirna"]="Overlap miRNAs";
//$queryTypes["pdomain"]="Overlap Protein Domains";
//$queryTypes["qtl"]="Overlap QTLs";
$queryTypes["genomeSNPValidator"]="Run genomeSNPValidator";
$queryTypes["genomeSNPAnnotator"]="Run genomeSNPAnnotator";
$queryTypes["getFastaSeq"]="Run getFastaSeq";
$queryTypes["getTxFastaSeq"]="Run getTxFastaSeq";
$functionClass["Exon Coding Synonymous"]="Coding Syn";
$functionClass["Exon Coding Non Synonymous"]="Coding nonSyn";

$GENOTYPEMAP{"A"}{"H"}="<td class='ah'>";
$GENOTYPEMAP{"A"}{"M"}="<td class='am'>";
$GENOTYPEMAP{"A"}{"L"}="<td class='al'>";
$GENOTYPEMAP{"T"}{"H"}="<td class='th'>";
$GENOTYPEMAP{"T"}{"M"}="<td class='tm'>";
$GENOTYPEMAP{"T"}{"L"}="<td class='tl'>";
$GENOTYPEMAP{"C"}{"H"}="<td class='ch'>";
$GENOTYPEMAP{"C"}{"M"}="<td class='cm'>";
$GENOTYPEMAP{"C"}{"L"}="<td class='cl'>";
$GENOTYPEMAP{"G"}{"H"}="<td class='gh'>";
$GENOTYPEMAP{"G"}{"M"}="<td class='gm'>";
$GENOTYPEMAP{"G"}{"L"}="<td class='gl'>";
$GENOTYPEMAP{"H"}{"H"}="<td class='h'>";
$GENOTYPEMAP{"N"}{"H"}="<td class='n'>";
/**************************************************************
This function filter the user input for security
***************************************************************/
function filterQuery($query){
  $query=htmlspecialchars(stripslashes($query));
  $query=str_ireplace("script", "blocked", $query);
  $query=mysql_escape_string($query);

 return $query;
}
/*****************************************************
This function displays the description/usage of
the selected query/tool
******************************************************/
function getXmlUrl($term){
 global $xml_file;$data="<div style='padding-top:2em;padding-left:2em;'>";
 $xml = simplexml_load_file($xml_file);
 if($xml){
   $tools=$xml->children();$i=0; 
   while($i<count($tools)){
     $name=$tools[$i]->attributes()->title;$tag=$tools[$i]->attributes()->name;
    if($tag!=$term){++$i;continue;}
     $data.="<h2>$name</h2><br/>";
     $notes=$tools[$i]->children()->note;$j=0;
     while($j<count($notes)){
           $note_name=$notes[$j]->attributes()->name;
           $data.="<h4>$note_name</h4>";
            $lists=$notes[$j]->children();
           if($notes[$j]->attributes()->id !=null){
             $data.="<ul>";$k=0;while($k<count($lists)){$data.="<li>$lists[$k]</li>";++$k;}
             $data.="</ul>";
           }
           else{$k=0;while($k<count($lists)){$data.="<p>$lists[$k]</p>";++$k;} }
          ++$j;
      }
     ++$i;
   }
 }
 $data.="</div>";
 return $data; 
}
/**************************************************************
 This function returns multiple selections of a given group
  into a nice formatted list where items are separated by ':'
***************************************************************/
function getSelectedItemList($myquerystring,$fieldname){
  $newlist=""; $i=0;
 if(!empty($myquerystring)){
    foreach($myquerystring as $term){ 
      if(preg_match("/$fieldname/",$term)){
         $fielddata=split("=",$term);
         $field_id=trim($fielddata[1]);
         if(!empty($field_id)){
            if($i==0)$newlist="$field_id";
            else{$newlist.=":$field_id";}
          }
          ++$i;
      }  
    }
  }
  return $newlist;
}
/************** getDefaultGenePrediction ***********
 get the list of annotation sources for this organism version
 and set the default annotation for this organism
****************************************************/
function getDefaultGenePrediction($orgv_id){
 $query="select gene_prediction_name from gene_prediction_by_organism op, gene_prediction p ";
 $query.="where op.organism_version_id =$orgv_id and op.gene_prediction_id=p.gene_prediction_id ";
 $query.=" and gene_prediction_name not in ('all_est','all_mrna','estOrientInfo')";
 global $con;if(!$con)$con=getConnection();
 $result = mysql_query($query,$con);$gene_prediction="";
 if($result){
    $numRows =mysql_numrows($result);$predictions=array();
   if($numRows>0){
       while($row = mysql_fetch_assoc($result)){
           $gene_prediction=$row["gene_prediction_name"];
           $predictions{$gene_prediction}=1;
       }
       if($predictions{"mgiGene"}==1){$gene_prediction="mgiGene";}
       elseif($predictions{"knownGene"}==1){$gene_prediction="knownGene";}
       elseif($predictions{"ensGene"}==1){$gene_prediction="ensGene";}
       elseif($predictions{"refGene"}==1){$gene_prediction="refGene";}
       elseif($predictions{"xenoRefGene"}==1){$gene_prediction="xenoRefGene";}
    }
  }
 return $gene_prediction;
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

/**************************************************
 displays the main page overview 
**************************************************/
function displayMain(){
 $data='<div class="about"><h2>About this tool</h2>
<p>Graber Transcriptdb is an attempt to create an integrated transcript database of all the organisms and annotations found in UCSC,MGI,ENSEMBL,UNIPROT,and EBI databases. 
 The database generates unique transcript_id and exon_id for all redundant transcripts and exons from various annotations for a given organism version.We store the start cordinates using the zero-base standard and the end cordinates using the one-base standard.
But we display all the cordinates using the one-base.  
</p>
<h2>You can use the tool as followed</h2>
 <ul><list>You have genomic regions of interest or a gene and you would like to extract overlaping functional features</li>
     <list>You have a file containing genomic regions of interest or a list of genes or transcripts and you would like to extract overlaping functional features</li>
 </ul>
 <p> Database updates are ran monthly to synchronize both our data and external data sources.
</p><hr/><h2>Good Practices:</h2>';
$data.=getXmlUrl("Good Practices");
$data.='</div>';
 return $data;
}
/******************************************************
 This function displays the search section of the page
*******************************************************/
function getToolDisplay($organism_version,$annotations,$query,$query_type){
 global $queryTypes; if(empty($query))$query="2:105509053-105543317";
 $data='<div id="tool-main"><div class="tool" style="z-index=0;"><table><tr>';
 //<td rowspan="2" width="70">         
 $options='<select name="queryType" id="anotselect" onchange="updateToolDoc();populateOrganismList();">';
 foreach ($queryTypes as $key=> $value){
      if($key=="$query_type")$options.="<option  value='$key' selected>$value</option>";
      else $options.="<option value='$key'>$value </option>";
  }
  $options.='</select><br/><br/>'; 
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
  $data.='onclick=\'this.value="";\'/><br/><br/>
          <h2>Or</h2><label for="file"><h3>upload a file Containing genomic cordinates:</h3></label><br/>
           <input type="file" name="file" id="file">
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
 $data.='<input type="submit" class="sb" name="submit" value="submit" onclick="validate();getFile();">
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
   global $txdb_service;$file=$txdb_service."?s=1&querytype=$querytype";$data="";
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
/********************************************************
 This function returns the list ofgene predictions
 for the specified organism version

*********************************************************/
function displayAnnotList($organism_version,$annotid){
   global $txdb_service;$file=$txdb_service."?v=$organism_version&l=1";
   $xml = simplexml_load_file("$file");$annot_list=array();$annot_map=array();
   $annot_list=explode(",",$annotid);foreach($annot_list as $annot)$annot_map["$annot"]=1;
   $data="<table><caption><div class='sh'><nobr>Available gene predictions for the selected organism version ($organism_version)</nobr></div></caption>";
   $data.="<tr>";
   if($xml){ $predictions=$xml->children()->prediction;$i=0;
      $line="";
       while($i<count($predictions)){  
          $annot_name=$predictions[$i]->children()->name;$this_annotid=$predictions[$i]->children()->gene_prediction_id;
          ++$i;
          if($annot_map["$annot_name"]==1)
             $line.="<li><input type='checkBox' name='annot[]' value=\"$annot_name\" checked/>&nbsp;$annot_name&nbsp;</li>";
          else $line.="<li><input type='checkBox' name='annot[]' value=\"$annot_name\"/>&nbsp;$annot_name&nbsp;</li>"; 
         if($i%4==0){
          $data.="<td><ul class='annotd'>$line</ul></td>"; $line="";
         }
      }
      if(!empty($line))$data.="<td><ul class='annotd'>$line</ul></td>";
      $data.="</tr></table>";
   }
   return $data;
}
function getFunctionalSNPsCount($query,$f_class,&$count,$source_id){
  global $snp_service;
   $url="$snp_service?format=xml&t=$query&fclass=$f_class&source=$source_id";global $functionClass;
   $xml = simplexml_load_file($url,null, LIBXML_NOCDATA);
   if($xml){
     $snps=$xml->children()->snp; $count=count($snps);
  }
}
/********************************************************
 This function returns the list of current organisms 
 that have gene annotations in our database
*********************************************************/
function getSNPsListing($query,$f_class,$line,&$count,$source_id,&$strainIndex){
   global $snp_service;
  /*Note: I will be getting SNPs from Ensembl snp web service :
     http://beta.rest.ensembl.org/vep/human/1:6524705:6524705/T/consequences?content-type=text/xml
     http://beta.rest.ensembl.org/documentation/user_guide
    using the chrom:chrStart:chrEnd:strand
  **/
   $url="$snp_service?format=xml&t=$query&fclass=$f_class&source=$source_id";global $functionClass;
  $xml = simplexml_load_file($url,null, LIBXML_NOCDATA);
  $rows="";
   if($xml){
     $snps=$xml->children()->snp; $index=0; $count=count($snps);
     $strains= getSNPsStrains($source_id);
     while($index < $count){$bg="";
         if($index%2==0)$bg="style='background:#eee;'";
         $snp=$snps[$index]; ++$index;
          $id= $snp->children()->id;$chrom= $snp->children()->chrom;
          $chrom.=":".$snp->children()->position;
          $alleles=$snp->children()->alleles;
          $genotype=$snp->children()->genotype;
          $mutation=$snp->children()->mutation; $ft_class=$snp->children()->function_class;
          $refcodon=$snp->children()->reference_aminoacid;$snpcodon=$snp->children()->snp_aminoacid;
          $snpOncds=$snp->children()->snp_posInCDS_posInProtein;
          $transcript=$snp->children()->transcript;$exonrank=$snp->children()->feature_rank;
          $gene=$snp->children()->gene;$ft_class=$functionClass["$ft_class"];
          $rows.="<tr $bg><td>$id</td><td>$chrom</td><td>$alleles</td>";
          $genotype_list=explode(";",$genotype);$strain_map=array();//$strains["$name"]["id"]=$id;
          foreach($genotype_list as $geno){
             list($allele,$strain_list)=explode(":",$geno); $strain=explode(",",$strain_list);
             foreach($strain as $str){ $id=$strains["$str"]["id"];
                 $strain_map{"$id"}=$allele;
              }
          }$geno="";global $GENOTYPEMAP;
          if(is_array($strainIndex)){$j=4;$conf="H";$char="N";
             while($j<count($strainIndex)){$id=$strainIndex[$j];++$j;
                $allele=$strain_map{"$id"};
                $geno.=$GENOTYPEMAP{"$allele"}{"$conf"}."$allele</td>";
             }
          }
           //<td>$genotype</td> ger
          $rows.="$geno<td>$mutation</td><td>$ft_class</td><td>$refcodon</td>";
          $rows.="<td>$snpcodon</td><td>$snpOncds</td><td>$exonrank</td><td>$transcript</td><td>$gene</td></tr>";
         
      }
   }
  return $rows; 
}
/********************************************************
 This function returns the list of current SNPs sources
 in our database
*********************************************************/
function getSNPsSources($source_id){
   global $snp_service;$url="$snp_service?format=xml&src=1";
   $xml = simplexml_load_file($url);$rows="<div class='sh'>Current SNPs Sources</div><ul id='sour'>";
   if($xml){
     $sources=$xml->children()->source; $index=0; $count=count($sources);
     while($index < $count){$bg="";
         if($index%2==0)$bg="style='background:#fff;'";
         $source=$sources[$index]; ++$index;$id= $source->children()->id;$name= $source->children()->name;
         if($id==$source_id) 
           $rows.="<li><input type='Radio' name='snpsource' value=\"$id\" checked/>&nbsp;$name&nbsp;</li>";
         else $rows.="<li><input type='Radio' name='snpsource' value=\"$id\"/>&nbsp;$name&nbsp;</li>";
         
      }
   }
  $rows.="</ul>";
  return $rows; 
}
/********************************************************
 This function returns the list of current strains
 and grayed strains that are not from the selected source
 
*********************************************************/
function getSNPsStrains($source_id){
   global $snp_service;$url="$snp_service?format=xml&str=1";
   $xml = simplexml_load_file($url,null, LIBXML_NOCDATA);
   $strains=array();
   if($xml){
     $sources=$xml->children()->strain; $index=0; $count=count($sources);
     while($index < $count){$bg="";
         $source=$sources[$index]; ++$index;$id= $source->children()->id;
         $name= $source->children()->name;
         $this_source=$source->children()->source;
         $strains["$name"]["id"]=$id;
         if($this_source==$source_id) $strains["$name"]["display"]=1;
         else $strains["$name"]["display"]=0; 
      }
  } 
  return $strains; 
}
/********************************************************
 This function displays the list of current strains
 and grayed strains that are not from the selected source
 
*********************************************************/
function displaySNPsStrains($source_id,$is_table,&$strainIndex){
     $strains=array(); $strains=getSNPsStrains($source_id);$rows="";$block="";$count=1;
     if($is_table!=1)$rows="<div class='sh'>Strains Selection</div><div id='str-list'>";
     $strain_keys=array_keys($strains); sort($strain_keys);$j=4;
     foreach($strain_keys as $name){
         $values=$strains["$name"]; $id=$values["id"]; 
         if($is_table!=1){ 
             $li="<li><input type='checkBox' name='strain' value=\"t$id\" id=\"$id\" 
                  onclick=\"toggleColumn(this,'snptable',1);\" disabled='true'/>&nbsp;$name&nbsp;</li>";
            if($values["display"]==1)
               $li="<li><input type='checkBox' name='strain' value=\"t$id\" id=\"$id\"  
                     onclick=\"toggleColumn(this,'snptable',1);\" />&nbsp;$name&nbsp;</li>";
            if($count%25==0){$rows.="<ul>$block</ul>";$block="";}
            else $block.=$li; ++$count;
         }else{
            $len=strlen($name);$verticalName=""; //<ul>";
            for($p=0;$p<=$len-1;++$p){
                if(!empty($name[$p]))$verticalName.="$name[$p]<br/>"; //</li>";
            }
            //$verticalName.="</ul>";
            $strainIndex[$j]=$id; $rows.="<th class='strainlabel' id='t$j'>$verticalName</th>";
           ++$j;
          
         }
      }
  if($is_table!=1){if(!empty($block))$rows.="<ul>$block</ul>";$block="";}
  if($is_table!=1)$rows.="</div>";
  return $rows; 
}
/********************************************************
 This function returns the list of transcripts that ovelap 
 the query
*********************************************************/
function getJasonTxListing($query,$organism_version,$annotations){
  global $txdb_service;$url="$txdb_service?t=$query&v=$organism_version&isoform=1";
  if(!empty($annotations))$url.="&a=$annotations";
 // echo "Calling $url\n";
  $exon_map=array();$tx_map=array();
  $xml = simplexml_load_file($url);$rows="";
   if($xml){
     $items= $xml->children()->annotation; $index=0; $count=count($items); 
     $head= $xml->children()->head;
     while($index < $count){
          $item=$items[$index]; ++$index;
          $chrom= $item->children()->chrom;$tx_start=$item->children()->chromStart;
          $tx_end=$item->children()->chromEnd;
          $chrom.=": $tx_start - $tx_end";$tx_len=$tx_end-$tx_start;
          $exon_len=0;
          $strand=$item->children()->strand;
          $transcript=$item->children()->transcriptNames;
          $gene=$item->children()->geneNames; $cds_start=$item->children()->cdsStarts;
          $cds_end=$item->children()->cdsEnds;
          $ex_starts=$item->children()->exonStarts;$ex_end=$item->children()->exonEnds;
     }
   }
}
/********************************************************
 This function returns the list of transcripts that ovelap 
 the query
*********************************************************/
function getTxListing($query,$organism_version,$isoform,&$count,&$head,$query_type,$annotations){
   global $txdb_service; 
   $url="$txdb_service?t=$query";if(!empty($organism_version))$url.="&v=$organism_version";
  if(!empty($isoform))$url.="&isoform=1";
  if(!empty($annotations))$url.="&a=$annotations";
 // echo "Calling $url\n";
  $exon_map=array();$tx_map=array();
  $xml = simplexml_load_file($url);$rows="";
   if($xml){
     $items= $xml->children()->annotation; $index=0; $count=count($items); 
     $head= $xml->children()->head;
     while($index < $count){$bg="";
         if($index%2==0)$bg="style='background:#eee;'";
          $item=$items[$index]; ++$index;
          $chrom= $item->children()->chrom;$tx_start=$item->children()->chromStart;
          $tx_end=$item->children()->chromEnd;
          $chrom.=": $tx_start - $tx_end";$tx_len=$tx_end-$tx_start;
          $exon_len=0;
          $strand=$item->children()->strand;
          $transcript=$item->children()->transcriptNames;
          $gene=$item->children()->geneNames; $cds_start=$item->children()->cdsStarts;
          $cds_end=$item->children()->cdsEnds;
          $ex_starts=$item->children()->exonStarts;$ex_end=$item->children()->exonEnds;
          $exon_count=explode(",",$ex_starts);$exonends=explode(",",$ex_end);$i=0;
          $tx_map{"$tx_start"}{"$tx_end"}=1;
          while($i<count($exon_count)){
            $exon_map{$exon_count[$i]}{$exonends[$i]}=1;++$i;
            $exon_len+=($exonends[$i]-$exon_count[$i]);
          }
          //for each cds start-end, find all the exons within the coding region then 
          // compute the len
          $cdsLen=getCDSlen($cds_start,$cds_end,$ex_starts,$ex_end); 
          $percenExon=$tx_len>0?$percenExon=round(($exon_len/$tx_len)*100,2):0;$f_class="1,2,5,6,7,8";$source_id=21;
          $percenCds=$exon_len>0?round(($cdsLen/$exon_len)*100,2):0;$snpcount=0;
          getFunctionalSNPsCount($transcript,$f_class,&$snpcount,$source_id);
          $rows.="<tr $bg><td><nobr>$chrom :$strand </nobr></td><td><nobr>$transcript - $gene</nobr></td>";
          $rows.="<td>$cds_start</td><td>$cds_end</td><td>".count($exon_count)."</td><td>$snpcount</td>";
          if($query_type=="exon")$rows.="<td>$ex_starts</td><td>$ex_end</td>";
          else $rows.="<td>$tx_len</td><td>$exon_len</td><td>$cdsLen</td><td>$percenExon</td><td>$percenCds</td>";
          $rows.="</tr>";
       }
      if($query_type=="exon")$count=count($exon_map);
      else $count=count($tx_map);
    }
  return $rows; 
}
function getCDSlen($cds_start,$cds_end,$ex_starts,$ex_end){
   $cdsSs=array();$exon_start=array();$exonends=array();$cdsLen=0;
   $cdsSs=explode(",", $cds_start);$cdsEs=explode(",",$cds_end);
   $exon_starts=explode(",",$ex_starts);$exonends=explode(",",$ex_end);$j=0;$cdscount=count($cdsSs);
   $ecount=count($exon_starts);
   while($j< $cdscount){$k=0;
        if(($cdsSs[$j]<=0)||(empty($cdsSs[$j]))){ ++$j; continue;}
        while($k < $ecount){
              if(($cdsSs[$j]>=$exon_starts[$k])&&($cdsSs[$j]<=$exonends[$k])){//fist coding exon
                     $cdsLen+=($exonends[$k]-$cdsSs[$j]);
              }else if(($cdsEs[$j]>=$exon_starts[$k])&&($cdsEs[$j]<=$exonends[$k])){
                      $cdsLen+=($cdsEs[$j]- $exon_starts[$k]);
             }else if(($exon_starts[$k]>$cdsSs[$j])&&($exonends[$k]<$cdsEs[$j])){
                $cdsLen+=($exonends[$k]-$exon_starts[$k]);} 
                 ++$k;  
          }
       ++$j;
   }
 return $cdsLen;
}
?>
