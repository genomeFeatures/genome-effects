<?php

/****************************************************************************
  Author: Lucie Ndzana Hutchins
  Date  : December 2012
  
  Note  : This utility for php function calls

*****************************************************************************/
require_once("user.php");

$dbname="cgd_snpdb";

/**************************************************************
This function filter the user input for security
***************************************************************/
function filterQuery($query){
  $query=htmlspecialchars(stripslashes($query));
  $query=str_ireplace("script", "blocked", $query);
  $query=mysql_escape_string($query);
  return $query;
}
/*********************************************************
 This function returns the chromosome id for
 a given chromsome name

**********************************************************/
function getChromosomeId($chrom){
 $query="select chromosome_id from snp_chromosome where chromosome_name='$chrom'";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$chr_id=0;
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
 This function returns the chromosome name for
 a given chromsome id

**********************************************************/
function getChromosomeName($chrom_id){
 $query="select chromosome_name from snp_chromosome where chromosome_id=$chrom_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$chr_name="";
 if($result){
    $numRows =mysql_numrows($result);
    if($numRows>0){
       $row = mysql_fetch_assoc($result);
       $chr_name=$row["chromosome_name"];
     }
  }
 return $chr_name;
}
/**********************************************************
 This function returns the a list of current SNPs sources
 found in our database
***********************************************************/
function getSourceList(){
 $query="select distinct s.source_id, source_name from snp_source s,";
 $query.=" snp_by_source sc where s.source_id=sc.source_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $sourceList=array();$source=array();
 $result = mysql_query($query,$cgd_con);$chr_name="";
 while($row= mysql_fetch_assoc($result)){
       $source[]=array('source'=> array('id'=>$row["source_id"],'name'=>$row["source_name"]));
  }
 $sourceList['sources']=$source;
 return $sourceList;
}
/**********************************************************
 this function returns the local id of a given source name
***********************************************************/
function getSourceId($soure_name){
 $query="select source_id from snp_source where source_name='$source_name'";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $id=0;$result = mysql_query($query,$cgd_con);
 if($result){$row= mysql_fetch_assoc($result);$id=$row["source_id"];}
 return $id;
}
/**********************************************************
 this function returns the name of a given strain id
***********************************************************/
function getStrainName($id){
 $query="select strain_name from snp_strain where strain_id=$id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $strain_name=0;$result = mysql_query($query,$cgd_con);
 if($result){$row= mysql_fetch_assoc($result);$strain_name=$row["strain_name"];}
 return $strain_name;
}
/**********************************************************
 This function returns the a list of strains for a given
 SNPs source
***********************************************************/
function getSourceStrainList($source_id){
 $query="select distinct s.strain_id,s.source_id,strain_name from snp_strain_by_source s,snp_strain sc ";
 $query.=" where s.strain_id=sc.strain_id ";
 if(!empty($source_id))$query.=" and s.source_id in($source_id)"; 
 $query.=" order by  s.strain_id";
 
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $strainList=array();$strain=array();
 $result = mysql_query($query,$cgd_con);$chr_name="";
 if(!$result)queryError($query,mysql_error()); //echo "$query";
 while($row= mysql_fetch_assoc($result)){ $name="<![CDATA[".$row["strain_name"]."]]>";
      $strain[]=array('strain'=> array('id'=>$row["strain_id"],'name'=>$name,'source'=>$row["source_id"]));
  }$strainList['strains']=$strain;
 return $strainList;
}
/**********************************************************************
 This function returns all the transcript that overlap a given genomic region
***********************************************************************/
function getTranscriptOverlap($chrom_id,$chrom_start,$chrom_end,$strand){
 $query="select t.transcript_local_id,t.gene_id,td.accession_id,chromosome_id,transcript_start,transcript_end,";
 $query.="transcript_strand from cgd_transcripts t,cgd_transcripts_desc td ";
 $query.=" where  td.transcript_local_id=t.transcript_local_id";
 $query.=" and chromosome_id=$chrom_id ";
 if(!empty($strand))$query.=" and transcript_strand='$strand'";
 $query.=" and transcript_end >= $chrom_start and transcript_start<=$chrom_end order by transcript_start";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 // if(!$result)queryError($query,mysql_error());
 while($row= mysql_fetch_assoc($result)){
      $gene_id=$row["gene_id"];
      $mgisymbol=getMGIGeneName($gene_id);if(!empty($mgisymbol))$mgisymbol=" ($mgisymbol)";
      $gene=getEnsGeneName($gene_id).$mgisymbol;
      $tx[]=array('transcript'=> array('id'=>$row["transcript_local_id"],'name'=>$row["accession_id"],
                   'gene'=>$gene,
                   'chromosome_id'=>$row["chromosome_id"],'start'=>$row["transcript_start"],
                   'end'=>$row["transcript_end"],'strand'=>$row["transcript_strand"]));
  }
 $txList['transcripts']=$tx;
 return $txList;
}

/**********************************************************************
 This function returns all the transcript that overlap a given genomic region
***********************************************************************/
function getTranscriptIDsOverlap($chrom_id,$chrom_start,$chrom_end,$strand){
 $query="select distinct transcript_local_id from cgd_transcripts ";
 $query.=" where chromosome_id=$chrom_id ";
 if(!empty($strand))$query.=" and transcript_strand='$strand'";
 $query.=" and transcript_end >= $chrom_start and transcript_start<=$chrom_end order by transcript_start";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 if(!$result)queryError($query,mysql_error());
 while($row= mysql_fetch_assoc($result)){
      $tx[]=$row["transcript_local_id"];
  }
 return $tx;
}
/**********************************************************************
 This function returns all the SNPs located on a give transcript/gene
***********************************************************************/
function getTxSNPs($tx_id,$gene_id,$function_class,$source_id,&$snpList,$geno){
 //$function_class="1,0,2";
 $query="select distinct t.snpid,t.transcript_local_id,t.gene_id,description,
         loc_rank from snp_transcript t, snp_loc_func l";
 $query.="  where t.snpid>0 ";
   if($tx_id>0)$query.=" and transcript_local_id=$tx_id ";
   if($gene_id>0)$query.=" and gene_id=$gene_id ";
 if(!empty($function_class)&&($function_class!=""))$query.=" and t._loc_func_key in($function_class)";
 $query.=" and t._loc_func_key=l._loc_func_key";
 if($tx_id<=0 && $gene_id<=0 && $chrom_start<=0)$query.=" limit 1000"; //default to the first 1000 SNPs
 $snps=array();$tx=array();//$snpList=array();
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 if(!$result)queryError($query,mysql_error());$table="snp_imputed";
 if(mysql_numrows($result)>0){
 while($row= mysql_fetch_assoc($result)){
       $mgisymbol=getMGIGeneName($row["gene_id"]);if(!empty($mgisymbol))$mgisymbol=" ($mgisymbol)";
       $gene=getEnsGeneName($row["gene_id"]).$mgisymbol;
       $snpid=$row["snpid"];$tx_id=getTranscriptTerm($row["transcript_local_id"]);
       $tx["$snpid"]["$tx_id"]["loc"]=$row["description"];
       $tx["$snpid"]["$tx_id"]["rank"]=$row["loc_rank"];
       $tx["$snpid"]["$tx_id"]["gene"]=$gene; $tx["$snpid"]["$tx_id"]["tx_id"]=$row["transcript_local_id"];
  }$snps_list=implode(",", array_keys($tx));$snpData=array();
  //print "got this list $snps_list for this transcript $tx_id ==\n";
  $snpData=getSNPsData($snps_list);$snps=getSNPsCordinates($snps_list);
  foreach ($snps as $snpid => $posData){
      $chrom_id=getChromosomeName($posData{"chrom"});$position=$posData{"pos"};
      $alleles=$snpData{"$snpid"}{"ref_allele"}."/".$snpData{"$snpid"}{"snp_allele"};$genotypes="";
      $mutation=$snpData{"$snpid"}{"mutation"};$snp_accessions=getSNPaccessionid($snpid,$source_id);
      if($geno){
         $genotype_array= getSNPsGenotype($snpid,$source_id,$table);
         foreach($genotype_array as $geno_allele => $strains){
           if($genotypes=="")$genotypes="$geno_allele:".implode(",",array_keys($strains));
           else $genotypes.=";$geno_allele:".implode(",",array_keys($strains));
         }
         $genotypes="<![CDATA[".$genotypes."]]>";
      }
      foreach ($tx["$snpid"] as $tx_id => $data){
           $snp_loc=$data{"loc"};$loc_rank=preg_match("/Intronic/i",$snp_loc)? "Intron ".$data{"rank"} :"Exon ".$data{"rank"};
           $transcript_id=$data{"tx_id"};
           $aa_implication=getSNPsImplication($snpid,$transcript_id);
           $pos="";$ref="";$other="";
           if(count($aa_implication)>0){
             $pos=$aa_implication["pos"];$ref=$aa_implication["ref"];$other=$aa_implication["other"];
           }   
           $snpList[]=array('snp'=>array('id'=>$snp_accessions,'chrom'=>$chrom_id,'position'=>$position,
             'alleles'=>$alleles,'genotype'=>$genotypes,'mutation'=>$mutation,'transcript'=>$tx_id,
              'function_class'=>$snp_loc,'feature_rank'=>$loc_rank,'snp_posInCDS_posInProtein'=>$pos,
              'reference_aminoacid'=>$ref,'snp_aminoacid'=>$other,'gene'=> $data["gene"]));
      }
  }
 }
}
/************************************************************
 Return the aminoacid implication of the snp on this transcript
*************************************************************/
function getSNPsImplication($snpid,$tx_id){
 $query="select concat(PosInCDS,':',PosInProtein) as pos,concat(ref_codon,':',ref_aa) as ref,";
 $query.="concat(snp_codon,':',snp_aa) as other from snp_aminoacid where snpid=$snpid and transcript_local_id=$tx_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$snpData=array();
 if($result){
   $row= mysql_fetch_assoc($result);
   $snpData{"pos"}=$row["pos"];$snpData{"ref"}=$row["ref"];$snpData{"other"}=$row["other"];
 }
 return $snpData;
}
/**********************************************************************
 This function returns the genomic location of each specified SNP id
***********************************************************************/
function getSNPsCordinates($snplist){
  $query="select snpid,chromosome_id, bp_position from snp_position where snpid in ($snplist)";
  $query.=" order by chromosome_id, bp_position ";$snps=array();
  global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
  $result = mysql_query($query,$cgd_con);
  while($row= mysql_fetch_assoc($result)){$snpid=$row["snpid"];
       $snps["$snpid"]["pos"]=$row["bp_position"];
       $snps["$snpid"]["chrom"]=$row["chromosome_id"];
  }return $snps;
}
/**********************************************************************
 This function returns the SNP alleles of each specified SNP id
***********************************************************************/
function getSNPsData($snplist){
  $query="select snpid,ref_allele,snp_allele,mutation_type ";
  $query.=" from snp_main where snpid in ($snplist)";
  $snps=array();
  global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
  $result = mysql_query($query,$cgd_con);
  if(!$result)queryError($query,mysql_error());
  while($row= mysql_fetch_assoc($result)){$snpid=$row["snpid"];
       $mutation=$row["mutation_type"]==1 ? "Transition":"Transversion";
       $snps["$snpid"]["ref_allele"]=$row["ref_allele"]; 
       $snps["$snpid"]["snp_allele"]=$row["snp_allele"];
       $snps["$snpid"]["mutation"]=$mutation;
  }return $snps;
}
/**********************************************************************
 This function returns the genotype alleles of every strain
 for the specified SNP id
***********************************************************************/
function getSNPsGenotype($snpid,$source_id,$table){
  $query="select genotype_allele,strain_id";
  $query.=" from $table where snpid=$snpid and source_id in($source_id) ";
  $snps=array();
  global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
  $result = mysql_query($query,$cgd_con);
  if(!$result)queryError($query,mysql_error());
  while($row= mysql_fetch_assoc($result)){$allele=strtoupper($row["genotype_allele"]);
        $id=$row["strain_id"];
        $strain=getStrainName($id);$snps["$allele"]["$strain"]=1;
  }return $snps;
}

function getSNPaccessionid($snpid,$source_id){
 $query="select distinct accession_id from snp_accession where snpid=$snpid ";
 //if(!empty($source_id))$query.=" and source_id in($source_id)";
 $query.=" and snpid_id_type in (1,3) order by snpid_id_type limit 1";
 $list=array();
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 while($row= mysql_fetch_assoc($result)){$list[]=$row["accession_id"];}
 //if(count($list)<=0)
 return implode(",", array_values($list));
}
function getStrainGeno($snpid,$source_id,$strain_id){
}
/*******************************************************
This function returns all the strains that contain the specified
genotype for a given SNP and source id
********************************************************/
function getStrainWithAllele($snpid,$allele,$source_id){
}
/**********************************************************************
this function returns the ensembl id of a given transcript
***********************************************************************/
function getTranscriptTerm($tx_id){
 $query="select accession_id from cgd_transcripts_desc where transcript_local_id=$tx_id";
 $id=""; global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 if($result){
    $row= mysql_fetch_assoc($result);$id=$row["accession_id"];
  }
 return $id;
}
/**********************************************************************
this function returns the first id of a given transcript term
***********************************************************************/
function getTranscriptID($tx_term){
 $query="select transcript_local_id from cgd_transcripts_desc where accession_id='$tx_term'";
 $id=0; global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);
 if($result){
    $row= mysql_fetch_assoc($result);$id=$row["transcript_local_id"];
  }
 return $id;
}
/**********************************************************************
this function returns the id and the cordinates of a given transcript
***********************************************************************/
function getTranscriptCordinates($tx_id,$gene_id){
 $query="select t.transcript_local_id,accession_id,chromosome_id,transcript_start,transcript_end,";
 $query.="transcript_strand from cgd_transcripts t,cgd_transcripts_desc td ";
 $query.=" where  td.transcript_local_id=t.transcript_local_id";
 if($tx_id>0)$query.=" and td.transcript_local_id=$tx_id ";
 if($gene_id>0)$query.=" and t.gene_id=$gene_id ";
 if($tx_id<=0 && $gene_id<=0)$query.=" limit 100"; //default to the first 100 transcripts
 $txList=array();$tx=array();
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$chr_name="";
// if(!$result)queryError($query,mysql_error());
 while($row= mysql_fetch_assoc($result)){
      $mgisymbol=getMGIGeneName($gene_id);if(!empty($mgisymbol))$mgisymbol=" ($mgisymbol)";
      $gene=getEnsGeneName($gene_id).$mgisymbol;
      $tx[]=array('transcript'=> array('id'=>$row["transcript_local_id"],'name'=>$row["accession_id"],
                   'gene'=>$gene,
                   'chromosome_id'=>$row["chromosome_id"],'start'=>$row["transcript_start"],
                   'end'=>$row["transcript_end"],'strand'=>$row["transcript_strand"]));
  }
 $txList['transcripts']=$tx;
 return $txList;
}
/**********************************************************************
this function returns the first id of a given gene accession_id
***********************************************************************/
function getGeneID($tx_term){
 $query="select gene_id from cgd_genes_desc ";
 $query.=" where accession_id='$tx_term'";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$id=0;
 if($result){
    $row= mysql_fetch_assoc($result);$id=$row["gene_id"];
  }
 if($id<=0){ //check if this is MGI gene
    $query="select gene_id from cgd_genes_ensembl_mgi ";
    $query.=" where marker_symbol='$tx_term' or mgi_geneid='$tx_term'";
    $result = mysql_query($query,$cgd_con);
    if($result){
       $row= mysql_fetch_assoc($result);$id=$row["gene_id"];
    }
  }
  if($id<=0){// check entrez gene table
    $query="select gene_id from cgd_genes_ensembl_entrez ";
    $query.=" where marker_symbol='$tx_term' or entrez_geneid='$tx_term'";
    $result = mysql_query($query,$cgd_con);
    if($result){
       $row= mysql_fetch_assoc($result);$id=$row["gene_id"];
    }
  }
 return $id;
}
/**********************************************************************
this function returns the id and the cordinates of a given gene
***********************************************************************/
function getGeneCordinates($gene_id){
 $query="select g.gene_id,gd.accession_id,chromosome_id,gene_start,gene_end,";
 $query.="gene_strand, gene_type from cgd_genes g,cgd_genes_desc gd ";
 $query.=" where g.gene_id=$gene_id and g.gene_id=gd.gene_id";
 $txList=array();$tx=array();
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$chr_name="";
 if(!$result)queryError($query,mysql_error());
 while($row= mysql_fetch_assoc($result)){
       $gene_strand=$row["gene_strand"]==1 ? "+" : "-";
       $mgisymbol=getMGIGeneName($gene_id);if(!empty($mgisymbol))$mgisymbol=" ($mgisymbol)";
       $gene=$row["accession_id"].$mgisymbol;
       $tx[]=array('gene'=> array('id'=>$row["gene_id"],'name'=>$gene,
                   'chromosome_id'=>$row["chromosome_id"],'start'=>$row["gene_start"],
                   'end'=>$row["gene_end"],'strand'=>$gene_strand));
  }
 $txList['genes']=$tx;
 return $txList;
}
/**********************************************************************
this function returns a list of transcripts of this gene
***********************************************************************/
function getGeneTxList($gene_id){
 $query="select distinct transcript_local_id from cgd_transcripts where gene_id=$gene_id";
 $txList=array();global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$chr_name="";if(!$result)queryError($query,mysql_error());
 while($row= mysql_fetch_assoc($result)){$txList[]=$row["transcript_local_id"];}
 return $txList;
}
/**********************************************************************
this function returns Ensembl annotations of this gene
***********************************************************************/
function getEnsGeneName($gene_id){
 $query="select accession_id from cgd_genes_desc where gene_id=$gene_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$gene_name="";
 if($result){$row= mysql_fetch_assoc($result);
    $gene_name=$row["accession_id"];
  }
 return $gene_name;
}
/**********************************************************************
this function returns MGI annotations of this gene
***********************************************************************/
function getMGIGeneName($gene_id){
 $query="select marker_symbol,mgi_geneid from cgd_genes_ensembl_mgi where gene_id=$gene_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$gene_name="";
 if($result){$row= mysql_fetch_assoc($result);
    $gene_name=$row["mgi_geneid"].",".$row["marker_symbol"];
  }
 return $gene_name;
}
/**********************************************************************
this function returns MGI annotations of this gene
***********************************************************************/
function getEntrezGeneName($gene_id){
 $query="select marker_symbol,entrez_geneid from cgd_genes_ensembl_entrez where gene_id=$gene_id";
 global $cgd_con;if(!$cgd_con)$cgd_con=getCgdConnection();
 $result = mysql_query($query,$cgd_con);$gene_name="";
 if($result){$row= mysql_fetch_assoc($result);
    $gene_name=$row["entrez_geneid"].",".$row["marker_symbol"];
  }
 return $gene_name;
}
/**********************************************************
 function to convert an associative array to xml format
***********************************************************/
function array_to_xml($annotations, &$xml) {
  foreach($annotations as $key => $value) {
     if(is_array($value)) {
        if(!is_numeric($key)){$subnode = $xml->addChild("$key");array_to_xml($value, $subnode);}
        else{ array_to_xml($value, $xml);}
     }
     else {$xml->addChild("$key","$value");}
   }
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
