 <?php /*********************************************************** 
  This file stores 
  all Global functions to be used accross the site
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : December 2009
  Modified : July 2010
  
  ***************************************************************/ 
   $newpath= ini_get('include_path');
   $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils"; 
   $newpath .=":/xraid/Users/LucieH/utils";
  
  ini_set('include_path', $newpath);
  require_once("build37/global.php");

   global $CGDTYPE;
   if(empty($CGDTYPE))loadCgdType();
 
/***** filters the input for funcky characters **/
function filterInput($data){
  if(!empty($data)){
   $data = preg_replace('/,/', '', $data); # remove the commas in case number input with commas
   $data = preg_replace("/'/", '', $data); # remove the quote
   $data = preg_replace('/"/', '', $data); # remove the quote
   $data = preg_replace('/</', '', $data); 
   $data = preg_replace('/\+/', '', $data); 
   $data = preg_replace("/>/", '', $data); 
   $data = preg_replace('/\?/', '', $data); 
  }
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
/****************************************************
* getSNPsource:
* Returns the source_name(s) of            
* selected  source_id(s)            
*****************************************************/
function getSNPsource($SNPsource){
   $snpsource=""; global $SNPSOURCE;
   $sourcedata=""; $is_imputed=0;
   $sourcelist=split(":",$SNPsource);
   foreach($sourcelist as $source_id){
        if(!empty($source_id)){
          $snpsource .= $SNPSOURCE{$source_id};
          $snpsource .=";";
          if(($source_id==16)||($source_id==21))$is_imputed=1;
        }
   } 
   $sourcedata{"sourcelist"}=$snpsource;
   $sourcedata{"imputedflag"}=$is_imputed; 
 return $sourcedata;
}

/************************************************************** 
 This function returns the error type if bad query
  or the genomic region chr:bpPosition if good query
****************************************************************/
function ValidateData($data){
   $error=""; $datatype=""; global $RANGEMAX; global $RANGEMAXs;global $Mbp;
   $data=trim($data);$data=preg_replace("/%3A/",":",$data);
  if(empty($data)){
     $datatype{"chr"}="";$datatype{"startpos"}="";$datatype{"endpos"}="";
     $datatype{"geneid"}=""; $datatype{"transcriptid"}="";
     $datatype{"qtype"}=-1;
     $datatype{"error"}="The '-- $data\n Enter Search term ' box should not be empty. Please enter a valid term and try again. "; 
  }
  else{
   $pos =preg_match("/^\s*(1|2|3|4|5|6|7|8|9|10|11|12|13|14|15|16|17|18|19|y|x|m(t?))\s*:\s*\d+(\.\d+)?/i",$data); 
     // check the format of the data
   $datatype{"query_single"}=$data;
   if($pos!=0) {           //Query using genome location
       list($chr,$position)=split(":",$data);
       $chr=trim($chr);
       if(preg_match("/x/i",$chr))$chr=20; 
       else if(preg_match("/y/i",$chr))$chr=21;
       else if(preg_match("/M/i",$chr)|| preg_match("/Mt/i",$chr))$chr=22;
       $datatype{"geneid"}=""; $datatype{"transcriptid"}="";
       $datatype{"qtype"}=0; $pos2=preg_match("/-/",$position);         // check if it is a range query 
       trim($position);
       if($pos2==0){       // this is a single position query
          $startpos= formatPosition($position);
          $endpos = $startpos;
       }
       else{             // this is a range query
          list($startpos,$endpos)=split("-",$position);
          trim($startpos); trim($endpos);
          $startpos =formatPosition($startpos);
          $endpos = formatPosition($endpos);
       }
       if($startpos<=0 ||$endpos<=0){ // bad query 
                                     //bad data format end pos is neither in bp nor in Mbp
          $datatype{"chr"}=$chr;$datatype{"startpos"}=$startpos;$datatype{"endpos"}=$endpos;
          $datatype{"error"}="Bad data format :[<b>$data</b>]. Start pos:$startpos; 
                      end pos:$endpos;Please review your query and try again."; 
       }
       else{
          if($endpos<$startpos){ // swap values if end_coord<start_coord
             $tempstr=$startpos;$startpos=$endpos;$endpos=$tempstr;
          }
          $datatype{"chr"}=$chr;$datatype{"startpos"}=$startpos;$datatype{"endpos"}=$endpos;
          $datatype{"error"}=""; 
          if(($endpos-$startpos)>$RANGEMAX){
              $datambp=($endpos-$startpos)/$Mbp;
              $datatype{"error"}="Query: <b>$data </b>( $datambp Mbp ) -- Genomic range too big. The maximum is $RANGEMAXs."; 
          }
        }
   }
   else{  //query using cgd type terms, now get the item type
      $cgdtypes=getCgdItemsType($con,$data); global $CGDTYPE;
      if(empty($cgdtypes)){
          $datatype{"chr"}=0;$datatype{"startpos"}=0;$datatype{"endpos"}=0;
          $datatype{"error"}=" We do not have SNP data for your search term:[<b> $data </b>] .please verify your query and try again."; 
          $datatype{"qtype"}=-1;
       } 
       else{
         //now use the in memory mapping to get the item info
         if(empty($CGDTYPE))loadCgdType();
         //GET the item_type_name,table_name, key_name
         $item_type_id=$cgdtypes[1]; $item=$cgdtypes[0];
         $item_type_name=$CGDTYPE{$item_type_id}{"type_name"};
         $item_table_name=$CGDTYPE{$item_type_id}{"item_table_name"};
         $item_key_name=$CGDTYPE{$item_type_id}{"item_key_name"};
         // now get the genomic region of this query
         $coord=getChrPosition($con,$item_type_id,$item);
         if($item_type_id==11){
             $source_id=getSNPAccessionSource($item);
             $datatype{"qsource"}=$source_id;
         }
         $chr= $coord[0];$startpos=$coord[1];$endpos= $coord[2];
         $gene_id=$coord[3];$transcript_id= $coord[4];
         $datatype{"chr"}=$chr;$datatype{"startpos"}=$startpos;$datatype{"endpos"}=$endpos;
         $datatype{"error"}=""; $datatype{"geneid"}=$gene_id;
         $datatype{"transcriptid"}=$transcript_id;
         $datatype{"qtype"}=$item_type_id;
       }  
   }
  }
   return $datatype;
}
//this function reformat the genome location in database format (int)
function formatPosition($postion){
  $pos=-1; global $Mbp;
  $pos =preg_match("/^\s*\d+\s*$/",$postion);
  if($pos !=0){ //position is decimal
     if($postion<10000){
       $pos=$postion*$Mbp; //convert position into Mbp if smaller int
     }
     else{
       $pos=$postion;
     }
  }
  else{  //check if position is in Mbp
     $pos=preg_match("/^\s*\d+\.\d+\s*$/",$postion);
     if($pos!=0){
        $pos=$postion*$Mbp;
     }
     else{  //this data has a bad format
       $pos=-1;
     }
  }
  return $pos;
}
/**** getCgdItemsType *********************************************
** This function returns the cgd type of a given search term *******
** CGD has the following types:
   Ensembl gene id   ->  7
   Ensembl transc id ->  8
   MGI gene id       ->  9
   MGI gene symbol   ->  10
   Good  SNP accession id  ->  11
******************************************************************/
function getCgdItemsType($con,$searchTerm){
   $searchTerm =trim($searchTerm);
   global $con;
   if(!$con)$con=getConnection();
   $query="select item,items_type_id from cgd_items";
   $query .=" where  item like '$searchTerm%' limit 1";
   $result = mysql_query($query,$con);
   $numRows =mysql_numrows($result);
   $item_type="";
   if($numRows>0){
      $row=mysql_fetch_assoc($result);
       $item_type[0]=$row["item"];
       $item_type[1]=$row["items_type_id"];
    }
    return $item_type;
}
/**** getStrainGroups *********************************** This function returns the list of all the strains
 by group
**************************************************************/
function loadStrainGroups(){
   global $con;global $STRAINGROUPS;
   if(!$con)$con=getConnection();
  $query="select strain_id,group_id from snp_strain_by_group ";
   $result = mysql_query($query,$con);$numRows =mysql_numrows($result);
   while($row=mysql_fetch_assoc($result)){
       $strain_id=$row{"strain_id"};$group_id=$row{"group_id"};
       $strain_id=trim($strain_id);$group_id=trim($group_id);
       $STRAINGROUPS{$group_id}{$strain_id}=1;
    }
 
}
/**** getThisSourceStrains *********************************** This function returns the list of all the strains
 provided by this source
**************************************************************/
function getThisSourceStrains($source_id,$group_id){
   global $con;
   if(!$con)$con=getConnection();
   $query="select s.strain_id from snp_strain_by_source s ";
   if(($group_id==1)||($group_id==4))$query.=",snp_strain_by_group g ";
   $query .=" where  s.source_id=$source_id";
  if(($group_id==1)||($group_id==4))$query.=" and g.group_id=$group_id and s.strain_id=g.strain_id ";
   $result = mysql_query($query,$con);
   //if(!$result) echo"Bad $query mysql_error()";
   $numRows =mysql_numrows($result);
   $strain_list=""; $i=0;
   while($row=mysql_fetch_assoc($result)){
       $strain_id=$row["strain_id"];
       $strain_id=trim($strain_id);
       if($i==0)$strain_list="$strain_id";
       else $strain_list.=":$strain_id"; ++$i;
    }
    return $strain_list;
}
/******getChrPosition ***************************************************
 Given a query term , Ensembl gene, transcript, MGI gene, SNP accession,
 ... this function does the mapping between 
 the term and the genome location (chromosome and bp_position
 *************************************************************************/
function getChrPosition($con,$type,$annotation){
   global $con;
   if(!$con)$con=getConnection();
  $annotation=trim($annotation);
  if(($type==10)||($type==9)){   //query by mgi gene symbol or id
     $query="SELECT g.chromosome_id,g.gene_start as start_loc,g.gene_end as end_loc,g.gene_id,";
     $query.=" '' as transc_id FROM cgd_genes_ensembl_mgi m, cgd_genes g ";
     $query .=" WHERE (m.marker_symbol= '$annotation' or m.mgi_geneid='$annotation')";
     $query .=" AND m.gene_id=g.gene_id ";
   }
   elseif($type==7){ //query by ensembl gene_id
     $query="SELECT g.chromosome_id,g.gene_start as start_loc,g.gene_end as end_loc,g.gene_id,";
     $query .=" '' as transc_id FROM cgd_genes g, cgd_genes_desc gd ";
     $query .=" WHERE gd.accession_id='$annotation' AND gd.gene_id=g.gene_id";
   }
   elseif(($type>10)&&($type<17)){ //query by snp_id
     $query="SELECT m.chromosome_id,m.bp_position as start_loc,m.bp_position as end_loc,'' as gene_id, ";
     $query .=" '' as transc_id FROM snp_accession a, snp_position m ";
     $query .=" WHERE a.accession_id='$annotation' AND a.snpid=m.snpid";
   }
   elseif(($type==17)||($type==18)){   //query by entrez gene
     $query="SELECT m.chromosome_id,m.gene_start as start_loc,m.gene_end as end_loc,m.gene_id,";
     $query.=" '' as transc_id FROM  cgd_genes m, cgd_genes_ensembl_entrez g ";
     $query.=" WHERE (g.marker_symbol='$annotation' or g.entrez_geneid='$annotation')";
     $query .=" AND g.gene_id=m.gene_id ";
   }
   elseif($type==19){   //query by ensembl protein
     $query="SELECT g.chromosome_id,p.cds_start as start_loc,p.cds_end as end_loc,g.gene_id,";
     $query.=" p.transcript_local_id as transc_id FROM cgd_ensembl_protein p,cgd_transcripts t, cgd_genes g ";
     $query .=" WHERE p.protein_id='$annotation'";
     $query .=" and p.transcript_local_id=t.transcript_local_id";
     $query .=" AND t.gene_id=g.gene_id ";
   }
   else{   //query using transcript
     $query="SELECT t.chromosome_id,t.transcript_start as start_loc,t.transcript_end as end_loc,t.gene_id,";
     $query .=" t.transcript_local_id as transc_id FROM cgd_transcripts_desc td, cgd_transcripts t ";// ,cgd_genes g ";
     $query .=" WHERE td.accession_id='$annotation' ";
  //   $query .=" AND t.gene_id=g.gene_id";
     $query .=" AND t.transcript_local_id=td.transcript_local_id ";
    }
   $result = mysql_query($query,$con); // or die("Bad query:$query".mysql_error());
   $numRows = mysql_numrows($result);
   if($numRows>0){
      $table="";
      $table[0]=mysql_result($result,0,"chromosome_id");
      $table[1]=mysql_result($result,0,"start_loc");
      $table[2]=mysql_result($result,0,"end_loc");
      $table[3]=mysql_result($result,0,"gene_id");
      $table[4]=mysql_result($result,0,"transc_id");
   }
  return $table;   
}
function setSNPloc($snp_loc){
    global $SNPLOC;
    $locdata[0]=0; $snp_location;$is_intergenic=0;
     $location_name=$SNPLOC{$snp_loc};
    $i=1;$snp_location[0]=$snp_loc;
    if(($snp_loc==0)||empty($snp_loc)){$snp_location[0]=0;$location_name="Intronic";}
    else if($snp_loc==12){
        $snp_location[0]=12;
        $location_name=$SNPLOC{12};
    }
    else if($snp_loc==10){
      $snp_location[0]=1;$snp_location[1]=2;$snp_location[2]=3;
      $snp_location[3]=4;$snp_location[4]=5;$snp_location[5]=6;$snp_location[6]=7; $snp_location[7]=8;$i=8;
    }
    else if($snp_loc==2){
      $snp_location[0]=2;$snp_location[1]=5; $snp_location[2]=6;
      $snp_location[3]=8;$snp_location[4]=7;$i=5;
    }
    else if(($snp_loc==34)||($snp_loc==3)){
      $snp_location[0]=3;$snp_location[1]=4; 
      $i=2;
    }
    else if($snp_loc==-1)$is_intergenic=1;
   $locdata[0]=$i;
   $locdata[1]= $snp_location;
   $locdata[2]=$location_name;
   $locdata[3]=$is_intergenic;
   
   return $locdata;
}
/***************************************************************
   Given the strain_list
   This function generates the list of strains with their IDs
****************************************************************/
function getSelectedStrainList($selectStrain_list){
        $strains=array();  global $con;$numRows=0;
        $strainIds=split(":",$selectStrain_list);
        if(is_array($strainIds)){
           if(!empty($strainIds[0])){
              if(!$con)$con=getConnection();
              $query="select distinct strain_id,strain_name from snp_strain ";
              $query.=" where strain_id in ("; $i=0;$strainid="";
              foreach($strainIds as $strain_id){
                 if($i==0){
                    if(!empty($strain_id)&&($strain_id>0))$strainid="$strain_id";
                  }
                 else{
                    if($strain_id>0)$strainid .=",$strain_id";
                  }
                 ++$i;
              }
              $query.="$strainid)";$query.=" order by strain_name";
              $result = mysql_query($query,$con) ;//or die("Bad query $query ".mysql_error());
              $numRows = mysql_numrows($result);
              $i=0;
              while($i<$numRows){
                    $strainID= mysql_result($result,$i,"strain_id");
                    $strainName = mysql_result($result,$i,"strain_name");
                    $strains{"$strainID"}=$strainName;
                    $i++;
              }
          }
     }
    return $strains;

}
/***************************************************************
   Given a snpid,geneid, transcriptid,
   This function generates the list of geneid where
   the SNP is located with their SNP function class
****************************************************************/
function getSNPlocation($snpid,$gene_id,$transcript_id,$snp_location){
   $query="select distinct gene_id,_loc_func_key";
   $query.=" from snp_transcript where snpid=$snpid ";
   if($gene_id>0)$query .=" and gene_id=$gene_id ";
   if($transcript_id>0)$query .=" and transcript_local_id=$transcript_id";
   //check the location 
   if(!empty($snp_location)){
      if(($snp_location!=-1)&& ($snp_location!=12)){
          if($snp_location==10)$query .=" and _loc_func_key in(1,2,3,4,5,6,7,8)";
          else if($snp_location==2) {
              $query .=" and _loc_func_key in(2,5,6,7,8)";
           }
           else $query.=" and _loc_func_key =$snp_location ";
       }
   }
   global $con;
   if(!$con)$con=getConnection();
   $result = mysql_query($query,$con);// or die("Bad query $query ".mysql_error());
   $numRows = mysql_numrows($result);
   $locdata="";
   $i=0;
   while($i<$numRows){
         $gene_id= mysql_result($result,$i,"gene_id");
         $loc_key=  mysql_result($result,$i,"_loc_func_key");
         $locdata{$gene_id} .="$loc_key:";
        ++$i;
   }  
  return $locdata; 
}

/***************************************************************
   Given a snpid
   This function generates the list of geneid, transcript where
   the SNP is located with their SNP function class
****************************************************************/
function getSNPGenelocation($snpid){
   $query="select distinct gene_id, transcript_local_id,_loc_func_key";
   $query.=" from snp_transcript where snpid=$snpid ";
   global $con;
   if(!$con)$con=getConnection();
   $result = mysql_query($query,$con);// or die("Bad query ".mysql_error());
   $numRows = mysql_numrows($result);
   $locdata="";
   $i=0;
   while($i<$numRows){
         $gene_id= mysql_result($result,$i,"gene_id");
         $transcript_id= mysql_result($result,$i,"transcript_local_id");
         $loc_key=  mysql_result($result,$i,"_loc_func_key");
         $locdata{$gene_id}{$loc_key} .="$transcript_id:";
        ++$i;
   }
  return $locdata; 
}

/****************************************************************
 Given a genomic region snp list and a filter set, this function
 returns a new list and number of SNPs where at least one of the selected
 strains has different genotype
 
*****************************************************************/
function getDifferenceList($SNPlist,$selectedStrains,$selectedSources,$is_imputed){
    $diffcount=0;$unique_source=array();
    foreach($sourcelist as $source_id){$unique_source{$source_id}=1;}
    $selectedSources="";
    while(list($source_id,$token)=each($unique_source)){
       if(empty($selectedSources))$selectedSources="$source_id";
       else $selectedSources.=":$source_id";
    }$sourcelist=split(":",$selectedSources); 
    $sourcecount=count($sourcelist);
    $impgenotypedata="";$genotypeData="";$i=0;$diffList[0]=0; $newList="";
    $snpcount=$SNPlist[0];$snpids=$SNPlist[1];
    while($i< $snpcount){
         $snpid=trim($snpids[$i]);
         if($is_imputed==1)$impgenotypedata=getGenotypedata($snpid,$selectedStrains,$selectedSources,$is_imputed);
         else{$genotypeData =getGenotypedata($snpid,$selectedStrains,$selectedSources,0);}
         $same_geno=1;
         if(!empty($genotypeData))$same_geno=$genotypeData{"samegeno"};
         else $same_geno=$impgenotypedata{"samegeno"};
         if($same_geno==0){$newList[$diffcount]=$snpid;++$diffcount;}
         ++$i;
     }
   $diffList[0]=$diffcount;
   $diffList[1]=$newList;
   return $diffList;
}
/****************************************************************
 this function returns the result table of a given query
 $currenttotalSNP,$strainsList,$SNPlist,$rowcount,$page,$datatype,$selectedStrains,
                        $is_imputed,$selectedSources,$strainIndex,$snp_loc,$fields,$querystring,$filterlink,$type
*****************************************************************/
function displayResulttable($currenttotalSNP,$strainsList,$SNPlist,$rowcount,$page,$datatype,$selectedStrains,
                        $is_imputed,$selectedSources,$strainIndex,$snp_loc,$fields,$querystring,$filterlink,$type){   
     global $strainIndex;$newSnpList=$SNPlist; $oldCount=$currenttotalSNP;//$oldCount=$SNPlist[0];
     $sg=$fields{"sg"};
     if($fields{"sg"}==1){
        $diffGenotypeList=getDifferenceList($SNPlist,$selectedStrains,$selectedSources,$is_imputed);
        $currenttotalSNP=$diffGenotypeList[0];//$SNPlist=$diffGenotypeList;
      }
     $data=displayPages($currenttotalSNP,$filterlink,$rowcount,$page,$snp_loc,$fields);
     $data.=getHeaderOptions($currenttotalSNP,$page,$rowcount,$fields,$querystring,$oldCount);
     if($type==0){
        $data.="<table border=1 id='snptable' width='99%' style='clear:both;'><thead style='font-size:0.90em;'>";
        $data.=getTableHeader($currenttotalSNP,$page,$rowcount,$strainsList,$fields,$querystring);
        $data.="</thead><tbody>";$startindex=0; 
        $data.=displayRows($SNPlist,$rowcount,$page,$datatype,$selectedStrains,$is_imputed,
             $selectedSources,$strainIndex,$snp_loc,$fields);
        $data.='<script>  allowDragging("snptable")  </script>';
        $data.="</tbody></table>"; 
     }
     else{
       $data.=getTableHeader($currenttotalSNP,$page,$rowcount,$strainsList,$fields,$querystring);
     }
    return "$data";                       
}
/****************************************************************
 this function returns the result filter options
 row count, filter genotype, download opntions
*****************************************************************/
function getHeaderOptions($currenttotalSNP,$page,$rowcount,$fields,$querystring,$oldCount){
  $i=1;$next=$currenttotalSNP; global $strainIndex; global $SNP_SEARCH_CSV;
     if(($page==1)||($page==0))$i=1; 
     else{$i=($page-1)*$rowcount;++$i;}
     $type=3;
     if(($i+$rowcount)<=$currenttotalSNP)$next=$i+$rowcount;   
     $diff="Results $i-$next of $currenttotalSNP";
     if($currenttotalSNP!=$oldCount)$diff="Result Total:$oldCount SNPs -->  $i-$next of $currenttotalSNP variations";
     $hint="<div id='cpage' style='width:90%; margin:auto;height:3em;clear:both;'>
            <div style='float:left;padding:0.4em;'>$diff</div>";
     $hint.="<div style=\"float:left;margin:auto;padding-left:1em;\">
                 <table><tr><td><b>Max Rows/Page :</b>
                           <select name=\"rcount\" id=\"rcount\" onchange=\"updatePageResult($type,0);\">
                            <option value=\"100\" selected=\"selected\">100</option>";
     if($rowcount==200)$hint.='<option value="200" selected="selected">200</option>';
     else $hint.='<option value="200">200</option>';
     if($rowcount==300)$hint.='<option value="300" selected="selected">300</option>';
     else $hint.='<option value="300" >300</option>';
     if($rowcount==500)$hint.='<option value="500" selected="selected">500</option>';
     else $hint.='<option value="500">500</option>';
     if($rowcount==1000)$hint.='<option value="1000" selected="selected">1000</option>';
     else $hint.='<option value="1000">1000</option>';
     $hint.='</select></td>';
     $checked=""; $type=1;
     if($fields{"sg"}==1)$checked="checked";                              
     $hint.="<td style='color:navy;padding-left:0.5em;'>
              <input type='checkbox' name='sg' id ='sg' onclick='updatePageResult($type,0);' value='1' $checked/>
                 Remove SNPs with no variation for selected strains 
             </td></tr></table>
              </div>
            <div style='float:right;padding:0.5em;margin:0;'>
            <a id='dload' href='$SNP_SEARCH_CSV?$querystring'><< Download Result >></a></div> </div>";
    return $hint;
}
function getTableHeader($currenttotalSNP,$page,$rowcount,$strainsList,$fields,$querystring){
     global $SNP_SEARCH_CSV;$type=0;global $strainIndex;
     
     $data="<tr><td class='sectionName' id='t1'>chrom</td><td class='sectionName' id='t2'>Position</td>";
     //now display strains header
      $j=3;$strainheader="";
     if(is_array($strainsList)){
        while(list($strain_id,$strain_name)=each($strainsList)){
           $len=strlen($strain_name);$verticalName="<ul>";
            for($p=0;$p<=$len-1;++$p){
                if(!empty($strain_name[$p]))$verticalName.="<li>$strain_name[$p]</li>";
            }
            $verticalName.="</ul>";
            $strainheader.="<td valign='bottom' class='strainlabel' id='t$j'>
                              $verticalName</td>";
            $strainIndex[$j]=$strain_id;
           ++$j;
         }
      }
      $field_name="SNPID";
      $len=strlen($field_name);$verticalName="";
      for($p=0;$p<=$len-1;++$p)$verticalName.="$field_name[$p]<br/>";
      $strainheader.=" <td class='sectionName' id='t$j'>$field_name</td>";++$j;
      $field_name="Variation Type";
      $len=strlen($field_name);$verticalName="";
      for($p=0;$p<=$len-1;++$p)$verticalName.="$field_name[$p]<br/>";
      if($fields{"v_type"}==1){
         $strainheader.="<td class='sectionName' id='t$j'>$field_name</td>";++$j;
      }
      if($fields{"g_name"}==1){
         $deleterow.="<td><a href='#' name=\"t$j\" onclick=\"toggleColumn(this,'snptable',$type);\" class='delete'/> x</a></td>";
         $strainheader.="<td class='sectionName' id='t$j'>MGI Gene</td>";++$j;
      }
      if($fields{"f_class"}==1){
         //$deleterow.="<td><a href='#' name=\"t$j\" onclick=\"toggleColumn(this,'snptable',$type);\" class='delete'/> x</a></td>";
         $strainheader.="<td class='sectionName' id='t$j' >Function class</td>";++$j;
      }
      if($fields{"cpg"}==1){
         //$deleterow.="<td><a href='#' name=\"t$j\" onclick=\"toggleColumn(this,'snptable',$type);\" class='delete'/> x</a></td>";
         $strainheader.="<td class='sectionName' id='t$j' >CpG</td>";++$j;
      }
      $strainheader.="<td class='sectionName' id='t$j' >Alleles</td>";++$j;
      if($fields{"m_allele"}==1){
         $strainheader.="<td class='sectionName' id='t$j'>Minor Allele %</td>";++$j;
      }
     if($fields{"m_genotype"}==1){
        //$deleterow.="<td><a href='#' name=\"t$j\" onclick=\"toggleColumn(this,'snptable',$type);\" class='delete'/> x</a></td>";
        $strainheader.="<td class='sectionName' id='t$j' title='Based on Max strain count'>Missing genotype %</td>";
        ++$j;
     }
     $strainheader.="<tr>";
      // now display the delete buttons    
      for($i=1; $i<$j;++$i){
          $strainheader.="<td>
           <a href='#' name=\"t$i\" onclick=\"toggleColumn(this,'snptable',$type);\" class='delete'/> x</a> </td>";
      }
      $strainheader="$strainheader</tr>";
      $data.="$strainheader";
    return $data;
}
/**************************************
Returns the genotype of this strain
****************************************/
function getThisGenotype($thisStraingeno,$html){
  global $GENOTYPEMAP; $geno="";$has_conflict=0;$delim=",";if(!$html)$delim="-";
  if(is_array($thisStraingeno)){
      list($first_allele,$confScore)=split("-",$thisStraingeno[0]);$geno_list=array();
      foreach($thisStraingeno as $allele){//get all the valid genotypes of this strain
         list($genotype_allele,$confScore)=split("-",$allele);
         if((strtoupper($genotype_allele)=="A")||(strtoupper($genotype_allele)=="T")||
            (strtoupper($genotype_allele)=="C")||(strtoupper($genotype_allele)=="G")){
           // if($geno_list{$genotype_allele}<=$confScore){
           $geno_list{$genotype_allele}=$confScore;//}
         }
      }
     if(count($geno_list)==0){$geno ="$has_conflict:<td class='n'>$first_allele</td>";}
     else{$temp_allele="";$uniq_allele=array();
         if(count($geno_list)>1)$has_conflict=1;$conf="H";$alleles="";$p=0;
         while(list($genotype_allele,$confScore)= each($geno_list)){
            if($p>0){$alleles.=",$genotype_allele";$uniq_allele{strtoupper($genotype_allele)}=1;}
            else{$alleles="$genotype_allele";
                $temp_allele=$genotype_allele;$uniq_allele{strtoupper($genotype_allele)}=1;++$p;}
            if(($confScore==1))$conf="M";else if($confScore==0)$conf="L";
         }
        if((strtoupper($temp_allele)=="A")||(strtoupper($temp_allele)=="T")||
           (strtoupper($temp_allele)=="C")||(strtoupper($temp_allele)=="G")){
            $temp_allele=strtoupper($temp_allele);$alleles="";
            $geno ="$has_conflict:".$GENOTYPEMAP{"$temp_allele"}{"$conf"};
            if(!$html)$geno ="$has_conflict:";
            while(list($this_allele,$score)=each($uniq_allele)){
              if(empty($alleles))$alleles=$this_allele;
              else $alleles.="$delim$this_allele";
            }
            if($html)$geno .="$alleles</td>";
            else{$geno .="$alleles";}
        }else{
          if($html)$geno .="$has_conflict:<td class='n'>$temp_allele</td>";
           else{$geno .="$has_conflict:$temp_allele";}
        }
     }
  }
  return $geno;
}
/***************************************************************
   Given the user query, this fuctions generates rows of
   SNPs located within the specified genomic region
****************************************************************/
function displayRows($snplist,$rowcount,$page,$datatype,$selectedStrains,
                     $is_imputed,$selectedSources,$strainIndex,$snp_location,$fields){
  $snpcount=$snplist[0];$snpids=$snplist[1];
  global $MGI_GENEDETAIL;global $SNP_SEARCH_HOME;global $CGD_GBROWSE_37;global $SNPLOC; 
  $chr=$datatype{"chr"};$startpos=$datatype{"startpos"};$endpos=$datatype{"endpos"};
  $gene_id= $datatype{"geneid"};global $GENOTYPEMAP;
  $transcript_id= $datatype{"transcriptid"}; $item_type_id=$datatype{"qtype"};
  $i=1; $sourcelist=split(":",$selectedSources); $sourcecount=count($sourcelist);
  if(($page==1)||($page==0))$i=0;
  else{$i=($page-1)*$rowcount;++$i;}$j=$i+$rowcount; $datarows="";
  $selectedStrainlist=split(":",$selectedStrains);$html=1;
  while(($i<$snpcount)&&($i<$j)){ // process every snp returned from selection
         $snpid=trim($snpids[$i]);$snpdata=getsnpdata($con,$snpid);
         $all=0;$snpaccession_id=getSNPaccession($snpid,$selectedSources,$all);
         if(!empty($snpdata)){
            $chr=$snpdata{"chr"};$bp_position=$snpdata{"pos"};
            $black6_allele=$snpdata{"b6allele"};$snp_allele=$snpdata{"snpallele"};
            $is_CpG_site=$snpdata{"iscpg"};$is_intergenic=$snpdata{"istergenic"};
            $cpg="no"; if($is_CpG_site)$cpg="yes";$mutation="Transversion";$mutation_type=$snpdata{"mutation_type"};
            if($mutation_type==1)$mutation="Transition";$impgenotypedata=array();$genotypeData=array();
            $mgigene="N/A";$ensemt="N/A";$snplocation="Intergenic";$is_same=0;$thisgeno="";
            if($is_intergenic==0)$snplocation="";$k=0; $l=3;
             $conf="H";$char="N"; $hasSameGenotype=0;
            $strains=split(":",$selectedStrains);
            if(!is_array($strains)&&(empty($strains[0]))){++$i;continue;}
            if($is_imputed==1)$impgenotypedata=getGenotypedata($snpid,$selectedStrains,$selectedSources,$is_imputed);
            else{$genotypeData =getGenotypedata($snpid,$selectedStrains,$selectedSources,0);}
            $same_geno=1;$strain_count=count($selectedStrainlist);
            if(!empty($genotypeData))$same_geno=$genotypeData{"samegeno"};
            else $same_geno=$impgenotypedata{"samegeno"};
            $actualcount=0;$thisgeno="";
            while($k<count($selectedStrainlist)){
                   $strain_id=$strainIndex[$l]; ++$l; ++$k;$conf="H";$char="N";//get this strain geno
                 if(!empty($impgenotypedata{"$strain_id"})){ //get the imputed genotype
                       $thisStraingeno=split(":",$impgenotypedata{"$strain_id"}); 
                       list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                       $thisgeno.=$geno;
                 }
                else if(!empty($genotypeData{"$strain_id"})){ // for other sources
                       $thisStraingeno=split(":",$genotypeData{"$strain_id"});
                       list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                       $thisgeno.=$geno;
                }
               else{// this strain does not have a genotype data /strain not provided by this source
                   $conf="H";$char="N";
                  if($strain_id==7){
                    $thisgeno.= $GENOTYPEMAP{"$black6_allele"}{"$conf"};
                    $thisgeno .="$black6_allele</td>";
                  }
                  else{$thisgeno .="<td class='n'>N</td>";} 
              }
          }//end of genotype list
         // if(count($selectedStrainlist)<=0){$thisgeno="";$same_geno=0;}
       /* */
      //////////////////////////////////////////////////////
      // get the genotype row for this snp
      // $thisgenotype=getGenotypedata();
      // check id all the strains have the same genotype
      // get minor allele frequency :
      /////////////////////////////////////////////////////
      $impgenoSum="";$genoSumarry="";
      //print "This is $snpid,$black6_allele,$selectedSources,$is_imputed\n";
      if($sourcecount==1)
         $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);
      else{
          $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,0);
          if($is_imputed==1)$impgenoSum= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);
       }
       $minorallele=""; $genomap="";
       $genomap= parseMinorAllele($genoSumarry,$impgenoSum,$snp_allele,$black6_allele);
       //now process the genomap to compute the minor allele frequency range.
       $minallelef=computeMAF($genomap,$black6_allele,$snp_allele,$sourcelist);
       $minal  =  $minallelef{"minal"};$minmaf =  $minallelef{"minmaf"};
       $nmap   =  $minallelef{"nmap"}; $totalStrain = $minallelef{"totalStrain"};
       $minorallele="$minal - $minmaf "; 
       if($is_intergenic==1){
          $row="<tr>"; 
          if($fields{"sg"}==1){if($same_geno==1) $row="<tr style='display:none;'>"; }                            
          $row.="<td class='tdata'>Chr$chr</td>";
          $row.="<td class='tdata'><a href=\"$CGD_GBROWSE_37/?name=$newchr:$gbrowse_start..$gbrowse_end\" target=\"_blank\"
                   title='CGD Gbrowse view of SNP:$newchr:$bp_position +-200 bp'>$bp_position </a></td>";
          $row.=$thisgeno; 
          $row.="<td class='tdata'><a href=\"$SNP_SEARCH_HOME?snpid=$snpid&ss=$selectedStrains\"
                    title='CGD SNP detail page for $snpaccession_id'>$snpaccession_id</a></td>";
          if($fields{"v_type"}==1){$row.="<td class='tdata'>$mutation</td>";}
          if($fields{"g_name"}==1){$row.="<td class='tdata'>$mgigene</td>";}
          if($fields{"f_class"}==1){$snplocation="Intergenic";
             $row.="<td class='tdata'><a href=\"$SNP_SEARCH_HOME?snpid=$snpid&ss=$selectedStrains\"
                    title='CGD SNP detail page for $snpaccession_id'>$snplocation</a></td>";
          }
          if($fields{"cpg"}==1){$row.="<td class='tdata'>$cpg</td>";}
          $row.=" <td class='tdata'>$black6_allele/$snp_allele</td>";
          if($fields{"m_allele"}==1){$row.="<td class='tdata'>$minorallele</td>";}
          if($fields{"m_genotype"}==1){$row .=" <td class='tdata'>$nmap</td>";}
          $datarows.="$row</tr>";
        }else{ $snplocation=""; 
             $snplocdata=getSNPlocation($snpid,$gene_id,$transcript_id,$snp_location);
            while(list($gene_id,$transcriptdata)=each($snplocdata)){
               $transcriptdatal=split(":",$transcriptdata);
               $gbrowse_start=$bp_position -200; $gbrowse_end=$bp_position+200;
               $newchr="Chr$chr";  $row="<tr>"; 
               if($fields{"sg"}==1) 
               if($same_geno==1) $row="<tr style='display:none;'>";                             
               $row.="<td class='tdata'>Chr$chr</td>";
               $row.="<td class='tdata'><a href=\"$CGD_GBROWSE_37/?name=$newchr:$gbrowse_start..$gbrowse_end\" target=\"_blank\"
                       title='CGD Gbrowse view of SNP:$newchr:$bp_position +-200 bp'>$bp_position </a></td>";
               $row.=$thisgeno; 
               $loccount=count($transcriptdatal);
               $row.="<td class='tdata'><a href=\"$SNP_SEARCH_HOME?snpid=$snpid&ss=$selectedStrains\"
                         title='CGD SNP detail page for $snpaccession_id'>$snpaccession_id</a></td>";
               if($fields{"v_type"}==1){$row.="<td class='tdata'>$mutation</td>";}
               if($fields{"g_name"}==1){
                  $mgigene=getMgiGene($gene_id,0);
                  $mgigene_symbol=$mgigene[1];$mgigene_id=$mgigene[0];
                  if(empty($mgigene))$mgigene_symbol=getEnsGene($gene_id);
                  $row.="<td class='tdata'><a href=\"$MGI_GENEDETAIL/accession_report.cgi?id=$mgigene_id\"";
                  $row.=" title='MGI Gene detail page for $mgigene_symbol' target='_blank'>$mgigene_symbol</a></td>";
                }$o=0;
                while($o<($loccount-1)){
                        $trans_location=$SNPLOC{$transcriptdatal[$o]};
                        if($transcriptdatal[$o]==1)$trans_location="CSyn";
                        else if($transcriptdatal[$o]==2)$trans_location="CNSyn";  
                        else if($transcriptdatal[$o]==5)$trans_location="StopGained";   
                        else if($transcriptdatal[$o]==6)$trans_location="StopLost";
                        else if($transcriptdatal[$o]==7)$trans_location="StartGained"; 
                        else if($transcriptdatal[$o]==0)$trans_location="Intronic";
                        else if($transcriptdatal[$o]==3)$trans_location="3'UTR";
                        else if($transcriptdatal[$o]==4)$trans_location="5'UTR";
                        else if($transcriptdatal[$o]==-3)$trans_location="Exon(non protein-coding)";
                        else if($transcriptdatal[$o]==-2)$trans_location="Intron(non protein-coding)";
                        else if($transcriptdatal[$o]==8)$trans_location="StartLost";                                        
                        if($o==0)$snplocation ="$trans_location";
                        else{$snplocation .=",$trans_location";} 
                        ++$o;                            
                 }
                if($fields{"f_class"}==1){$row.="<td>$snplocation</td>";}
                if($fields{"cpg"}==1){$row.="<td class='tdata'>$cpg</td>";}
                $row.=" <td class='tdata'>$black6_allele/$snp_allele</td>";
                if($fields{"m_allele"}==1){$row.="<td class='tdata'>$minorallele</td>";}
                if($fields{"m_genotype"}==1){$row .=" <td class='tdata'>$nmap</td>";}
                $datarows.="$row</tr>";
            }//end of  while(list($gene_id,$transcriptdata)=each($snplocdata)){
          } //end of if intergenic==0
        // }// while($k<count($selectedStrainlist)){
     }//end of if(!empty($snpdata)){
    ++$i;
  }//end of while(($i<$snpcount)&&($i<$j)){
  return $datarows;
}
function getSNP_coordinates($snpid){
  global $con;if(!$con)$con=getConnection();
  $query="select chromosome_id,bp_position from snp_position where snpid=$snpid";
  $result = mysql_query($query,$con); // or die ("Bad query:$query".mysql_error());^M
  $numRows = mysql_numrows($result);$chrom="";
  if($numRows>0){
     $chrom{"pos"}=mysql_result($result,0,"bp_position");
     $chrom{"chr"}=mysql_result($result,0,"chromosome_id");
  }
  return $chrom;
}
/*********************************************************
 getsnpdata:
 Given a snpid,returns:
 chromosome_id,bp_position,snp_allele,black6_allele,
 is_conflict,is_intergenic
 **********************************************************/
function getsnpdata($con,$snpid){
    $snpdata="";$snpid=trim($snpid);
    global $con;if(!$con)$con=getConnection();
  if($snpid>0){
    $query="select snp_allele,ref_allele,is_intergenic,is_CpG_site,mutation_type,is_conflict";
    $query.=" from snp_main where snpid=$snpid";
    $result = mysql_query($query,$con); // or die ("Bad query:$query".mysql_error());
    $numRows = mysql_numrows($result);
    if($numRows>0){
        $chrom=getSNP_coordinates($snpid);
        $black6_allele = mysql_result($result,0,"ref_allele");
        $snp_allele= mysql_result($result,0,"snp_allele");
        $is_CpG_site= mysql_result($result,0,"is_CpG_site");
        $is_intergenic= mysql_result($result,0,"is_intergenic");
        $mutation_type= mysql_result($result,0,"mutation_type");
        $is_conflict=mysql_result($result,0,"is_conflict");
        $snpdata{"pos"}=$chrom{"pos"};
        $snpdata{"b6allele"}=$black6_allele;
        $snpdata{"snpallele"}=$snp_allele;
        $snpdata{"iscpg"}= $is_CpG_site;
        $snpdata{"chr"}= $chrom{"chr"};
        $snpdata{"istergenic"}= $is_intergenic;
        $snpdata{"mutation_type"}= $mutation_type;
        $snpdata{"is_conflict"}=$is_conflict;
    }
  }
    return $snpdata;
}
/*************************************************************************************
 Given a genomic region (chr, start, end) /and snp location (intron, exon, 3'utr,..)
 returns the list of snpids
**************************************************************************************
function getRangeSNPList($chr,$startpos,$endpos,$loc){
 $query="select distinct m.snpid from snp_main m ";
 $loc="($loc)";
 if($loc==10)$loc="(1,2,3,4,5,6,8)";
 else if($loc==2)$loc="(2,5,6,8)";
 if($loc!=12){
    $query.=" ,snp_transcript t";
    $query.=" where m.chromosome_id=$chr and m.is_intergenic=0  and m.bp_position between $startpos and $endpos";
    $query.=" and m.is_public=1 and m.snpid=t.snpid and t._loc_func_key   in $loc";
  }
  else{
    $query.=" where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos";
    $query.=" and m.is_public=1";
  }
  global $con;
  if(!$con)$con=getConnection(); $newlist="";$thisregionsnplist[0]=0;$thisregionsnplist[1]="";
  $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
  $numRows = mysql_numrows($result); $i=0;
  while($row=mysql_fetch_assoc($result)){
          $snpid=$row["snpid"]; 
          $newlist[$i]=$snpid; 
          
          ++$i;     
  }
  $thisregionsnplist[0]=$numRows;
  $thisregionsnplist[1]=$newlist;
  return $thisregionsnplist;
 
}
/*********************************************************************** 
 create a memory structure for SNPs within a given genomic region
// to include: source list, strains list, is intergenic flag,loc_list
************************************************************************/
function loadQuerySNPsInSession($chr,$startpos,$endpos){
  global $con; $allloc=12;$intronic=0;$exonall=10; $urt_3=3; $utr_5=4;
  $synonymous=1; $non_syn=2;$Stop_Gained=5;$Stop_Lost=6;$Initial_Met=8;
  if(!isset($_SESSION['genome']))$_SESSION['genome']=array();
  if(!isset($_SESSION['snp']))$_SESSION['snp']=array();
  $allSNPs=array();
  if(!$con)$con=getConnection(); $newlist="";
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$allloc"])){
     $allSNPs=getRangeSNPList($chr,$startpos,$endpos,$allloc);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$allloc"]["snpCount"]=$allSNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$allloc"]["snpList"]=$allSNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$intronic"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$intronic);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$intronic"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$intronic"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$exonall"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$exonall);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$exonall"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$exonall"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$urt_3"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$urt_3);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$urt_3"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$urt_3"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$utr_5"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$utr_5);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$utr_5"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$utr_5"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$synonymous"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$synonymous);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$synonymous"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$synonymous"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$non_syn"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$non_syn);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$non_syn"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$non_syn"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Gained"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$Stop_Gained);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Gained"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Gained"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Lost"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$Stop_Lost);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Lost"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Stop_Lost"]["snpList"]=$SNPs[1];
  }
  if(!isset($_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Initial_Met"])){
     $SNPs=getRangeSNPList($chr,$startpos,$endpos,$Initial_Met);
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Initial_Met"]["snpCount"]=$SNPs[0];
     $_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$Initial_Met"]["snpList"]=$SNPs[1];
  }
  // now for each SNP id from the main list, get supporting snp data
  foreach($allSNPs[1] as $snpid){
      if(!isset($_SESSION['snp']["$snpid"])){
         // get this snp data
         $snpdata=getsnpdata($con,$snpid);
         $_SESSION['snp']["$snpid"]['pos']=$snpdata{"pos"};
         $_SESSION['snp']["$snpid"]['b6allele']=$snpdata{"b6allele"};
         $_SESSION['snp']["$snpid"]['snpallele']=$snpdata{"snpallele"};
         $_SESSION['snp']["$snpid"]['iscpg']=$snpdata{"iscpg"};
         $_SESSION['snp']["$snpid"]['chr']=$snpdata{"chr"};
         $_SESSION['snp']["$snpid"]['istergenic']=$snpdata{"istergenic"};
         $_SESSION['snp']["$snpid"]['mutation_type']=$snpdata{"mutation_type"};
         $_SESSION['snp']["$snpid"]['is_conflict']=$snpdata{"is_conflict"};
         $_SESSION['snp']["$snpid"]['loc']=array();
         $_SESSION['snp']["$snpid"]['transcript']=array();
         $_SESSION['snp']["$snpid"]['gene']=array();
         //now get all the genomic locations for this snp
         if($snpdata{"istergenic"}==0){
            $_SESSION['snp']["$snpid"]['loc']= getThisSNPLocation($snpid);
            $_SESSION['snp']["$snpid"]['transcript']= getThisSNPTranscript($snpid);
            $_SESSION['snp']["$snpid"]['gene']= getThisSNPGene($snpid);
          }
         // now get all the sources of this snp
        $_SESSION['snp']["$snpid"]['source']= getThisSNPSource($snpid);
        
      }
  }
    
}
/**************************************************************************
 For get all the sources of this SNP
***************************************************************************/
function getThisSNPSource($snpid){
  $query="select distinct source_id from snp_by_source where snpid=$snpid";
  global $con; $sourcelist=array();$i=0;
  if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
  while($row=mysql_fetch_assoc($result)){
          $source_id=$row["source_id"];
          $sourcelist{"$source_id"}=$source_id;  
  }
  return $sourcelist;
}
/**************************************************************************
 For SNP located on transcript, get all the functional implications
  of this SNP
***************************************************************************/
function getThisSNPLocation($snpid){
  $query="select distinct _loc_func_key from snp_transcript where snpid=$snpid";
  global $con; $loclist=array();$i=0;
  if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
  while($row=mysql_fetch_assoc($result)){
          $loc_id=$row["_loc_func_key"];
          $loclist{"$loc_id"}=$loc_id;     
  }
  return $loclist;
}
/**************************************************************************
 For SNP located on transcript, get all the transcripts
  of this SNP
***************************************************************************/
function getThisSNPTranscript($snpid){
  $query="select distinct transcript_local_id from snp_transcript where snpid=$snpid";
  global $con; $loclist=array();$i=0;
  if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
  while($row=mysql_fetch_assoc($result)){
          $loclist[$i]=$row["transcript_local_id"]; 
          ++$i;     
  }
  return $loclist;
}
/**************************************************************************
 For SNP located on transcript, get all the transcripts
  of this SNP
***************************************************************************/
function getThisSNPGene($snpid){
  $query="select distinct gene_id from snp_transcript where snpid=$snpid";
  global $con; $loclist=array();$i=0;
  if(!$con)$con=getConnection();
  $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
  while($row=mysql_fetch_assoc($result)){
          $loclist[$i]=$row["gene_id"]; 
          ++$i;     
  }
  return $loclist;
}
/*************************************************************************** 
   getSNPsIDs: given a genomic region,source(s),and annotation(s)
    this function returns a list of SNPs Ids within the specified region
*****************************************************************************/
function getSNPsIDs($con,$chr,$startpos,$endpos,$selectedSources,$snp_location,$is_intergenic,
        $gene_id,$transcript_local_id,$locCount){
   global $con; $thisregionsnplist=""; global $SNPLIST;
   if(!$con)$con=getConnection(); $newlist="";
    // get the user selected sources
   $sources= split(":",$selectedSources); 
   $query="select distinct sb.snpid";
   $query .=" from snp_by_source sb , snp_position m ";
   if(($gene_id>0)||($transcript_local_id>0)){ //query using gene annotations
       $query.=",snp_transcript st where  m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
       $query.=" and st.snpid=m.snpid  ";
       if($gene_id>0) $query.=" and st.gene_id=$gene_id ";
       if($transcript_local_id>0)$query.=" and st.transcript_local_id=$transcript_local_id ";
       //check the location 
       if(($snp_location[0]!=-1)&&$snp_location[0]!=12){
             $i=0; $query .=" and st._loc_func_key in(";
             while($i<$locCount){
                   if($i==0)$query .="$snp_location[$i]";
                   else{
                       if($snp_location[$i]>0)$query .=",$snp_location[$i]";
                     }
                     ++$i;
              }
              $query .=") ";
        }
    }
    else{
       if($is_intergenic==1){
          $query.=" ,snp_main ms where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
          $query .=" and m.snpid=m.snpid and ms.is_intergenic=1 ";
        }
       else{
           if(($snp_location[0]!=-1)&&($snp_location[0]!=12)){
               $query.=",snp_transcript st ";
               $query.=" where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
               $query .=" and m.snpid=st.snpid ";
               $i=0; $query .=" and st._loc_func_key in(";
               while($i<$locCount){
                     if($i==0)$query .="$snp_location[$i]";
                     else{
                        if($snp_location[$i]>0)$query .=",$snp_location[$i]";
                     }
                     ++$i;
               }
               $query .=") ";
           }
           else{
               $query.=" where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
           }
         }
    }
    $query .=" and m.snpid=sb.snpid and sb.source_id in(";
    $i=0;
    foreach($sources as $source_id){
           if($i==0)$query .="$source_id";
           else{
              if($source_id>0)$query .=",$source_id";
            }
          ++$i;
    }
    $query .=") ";
    /*$tx_table="";$snp_main="";
    if($snp_location==-1)$snp_main=", snp_main ms";
    if(($gene_id>0)||($transcript_local_id>0)||
       (($snp_location!=-1)&&($snp_location!=12)))$tx_table=",snp_transcript ts";
    $query="select distinct p.snpid from snp_position p $tx_table, snp_by_source t";
    $query.=" where p.chromosome_id=$chr and bp_position between $startpos and $endpos";
    $query.=" and p.snpid=t.snpid and t.source_id in(";
    $sources= split(":",$selectedSources);$i=0;
    foreach($sources as $source_id){
       if($i==0)$query .="$source_id";
       else{if($source_id>0)$query .=",$source_id";}++$i;
    }
    $query.=")";$term="ts.transcript_local_id=$transcript_local_id";
    if(($gene_id>0)||($transcript_local_id>0)||
       (($snp_location!=-1)&&($snp_location!=12)))$query.=" and p.snpid=ts.snpid ";
    if($gene_id>0)$query.=" and ts.gene_id=$gene_id ";if($transcript_local_id>0)$query .=" and $term ";
    if(($snp_location!=-1)&& ($snp_location!=12)){
       if($snp_location==10){$query .=" and ts._loc_func_key in(1,2,3,4,5,6,7,8)";}
       else if($snp_location==2){$query .=" and ts._loc_func_key in(2,5,6,7,8)";}
       else if($snp_location==34){$query .=" and ts._loc_func_key in(3,4)";}
       else $query.=" and st._loc_func_key =$snp_location ";
    }*/
    $query.=" order by m.bp_position";
    $result = mysql_query($query,$con);// or die ("Bad query:$query".mysql_error());
    $numRows = mysql_numrows($result); $i=0;
    while($row=mysql_fetch_assoc($result)){
          $snpid=$row["snpid"]; 
          $newlist[$i]=$snpid; 
          ++$i;     
    }
    $thisregionsnplist[0]=$numRows;
    $thisregionsnplist[1]=$newlist;
    return $thisregionsnplist;
}
/****/
function getOtherGeno($snpid,$table_name,$strain_id){
 $query="select distinct genotype_allele from $table_name  where snpid=$snpid and strain_id=$strain_id and genotype_allele in ('A','T','C','G')";
 global $con;$temp_allele="";
 if(!$con)$con=getConnection();$result = mysql_query($query,$con);
 //if(!$result)echo "$query";
  while($row=mysql_fetch_assoc($result)){ $geno=$row["genotype_allele"];
     if(empty($temp_allele))$temp_allele="$geno";
     else $temp_allele.=",$geno";
  }
 return $temp_allele;
}
/*******************************************************************
  Given a snpid, list of strains, and a list of sources
  this function returns the genotype allele information 
  including the samegeno flag set to 0 if at least one of the strains
  has different genotype
********************************************************************/
function getGenotypedata($snpid,$selectedStrains,$selectedSources,$is_imputed){
   $strain_genotypes="";  $sources="";$numRows=0;$strains=split(":",$selectedStrains);
   $listOfSelectedStrains=array();
   if(!empty($selectedSources))$sources=split(":",$selectedSources);
   $strain_genotypes=array();
   if(is_array($strains)){//initiate all the selected strains with genotype=N
       foreach($strains as $strain_id){
           if($strain_id>0){
             $listOfSelectedStrains{"$strain_id"}=1;$strain_genotypes{"$strain_id"}="";}
       }
    }
   $query="select distinct strain_id,genotype_allele from snp_genotype ";
    if($is_imputed)
      $query="select strain_id,genotype_allele,confidence from snp_imputed";
   $query .=" where snpid=$snpid  ";
    if(!empty($selectedSources)){
         $query.=" and source_id in(";
         foreach($sources as $source_id){
            if($source_id>0){
               if($i==0){$query.="$source_id";}
               else{$query.=",$source_id";}
            }
            ++$i;
          }
        $query .=") ";
       }
   $query .=" order by strain_id ";
    global $con;  $i=0;$geno_list=array();$temp_allele="";
    if(!$con)$con=getConnection();$result = mysql_query($query,$con);
    if($result){$prevallele=""; $geno_list="";$strain_genotypes{"samegeno"}=1;$geno_list="";
        while($row=mysql_fetch_assoc($result)){
             $strainID= $row["strain_id"];$genotype_allele = $row["genotype_allele"];
             if(empty($listOfSelectedStrains{"$strainID"}))continue;
             if(!$is_imputed){
                  $other_table="snp_imputed";$other_genos=getOtherGeno($snpid,$other_table,$strainID);
                  if(empty($other_genos))
                     $strain_genotypes{"$strainID"}.="$genotype_allele-2:"; //high confidence score=2
                  else{ $genos=split(",",$other_genos); $otherallele="";
                     foreach($genos as $geno){
                        if(strtoupper($genotype_allele)!=strtoupper($geno)){$otherallele=$geno;}
                     }
                     if(empty($otherallele))
                        $strain_genotypes{"$strainID"}.="$genotype_allele-2:"; 
                     else $strain_genotypes{"$strainID"}.="$genotype_allele-2:$otherallele-2:";
                  }
             }
             else{ $confidence= $row["confidence"];
                  $other_table="snp_genotype";$other_genos=getOtherGeno($snpid,$other_table,$strainID);
                   if(empty($other_genos))
                      $strain_genotypes{"$strainID"}.="$genotype_allele-$confidence:";
                    else{ $genos=split(",",$other_genos); $otherallele="";
                          foreach($genos as $geno){
                              if(strtoupper($genotype_allele)!=strtoupper($geno)){$otherallele=$geno;}
                           }
                          if(empty($otherallele))$strain_genotypes{"$strainID"}.="$genotype_allele-$confidence:";
                           else {$geno_list{$otherallele}=1;
                              $strain_genotypes{"$strainID"}.="$genotype_allele-$confidence:$otherallele-$confidence:";
                           }
                    }
              }
             if((strtoupper($genotype_allele)=="A")||(strtoupper($genotype_allele)=="T")||
                 (strtoupper($genotype_allele)=="C")||(strtoupper($genotype_allele)=="G")){
                $temp_allele=strtoupper($genotype_allele);$geno_list{$temp_allele}=1;
             }
         }
         // now set the same genotype flag
         if(count($geno_list)>1)$strain_genotypes{"samegeno"}=0;
     }
     return  $strain_genotypes;
}

/****************************  getSNPaccession ********************
 Returns a list of all the accession ids of this snp
 ******************************************************************/
function getSNPaccession($snpid,$selectedSources,$all){
         $i=0;global $con;if(!$con)$con=getConnection();$accession="";
        $query =" select accession_id from snp_accession where snpid=$snpid and snpid_id_type not in(2,4)";
        if($all<=0){
          $sources=split(":",$selectedSources);
           $query.=" and source_id in(5";
           foreach($sources as $source_id){
                  if($source_id>0){
                     $query.=",$source_id";
                   }
            }
           $query .=") ";
        }
        $query.=" order by snpid_id_type ";
        if($all==0)$query .=" limit 1";
       // print"The query is :$query -----";
        $result = mysql_query($query,$con);// or die ("Bad query:$query");
        $numRows = mysql_numrows($result);
        if($numRows>0){
           $i=0;
           if($numRows==1)$accession = mysql_result($result,$i,"accession_id");
           else{
              while($i<$numRows){
                    $accession .= mysql_result($result,$i,"accession_id");
                    if($i<$numRows-1)$accession .=",";
                    ++$i;
              }
          }
        }
        return $accession;
}
/***************************************************************
  given an aminoacid, this function returns the frequency 
  of the the codon with max codon_fraction_1000 
  *************************************************************/
function getFmax($codon,$con){
  $aa= getAminoacid($codon,$con);
  $query="select max(codon_fraction_1000) as f_max from cgd_mouse_codonFrequency  ";
  $query .=" where aa_one_letter_abrev='$aa' and organism_id=1";
  $result = mysql_query($query,$con);// or die("Bad query:$query");
  $numRows =mysql_numrows($result);
  $fmax=0;;
  if($numRows>0){
     $row=mysql_fetch_assoc($result);
     $fmax=$row["f_max"];
  }
  return $fmax;
}
/***************************************************************
  given a codon, this function returns the corresponding 
  one letter abreviation of the aminoacid of the the codon 
  *************************************************************/
function getAminoacid($codon,$con){
  $newcodon=dna2rna($codon);
  $query="select aa_one_letter_abrev as aa from cgd_mouse_codonFrequency";
  $query .=" where codon='$newcodon' and organism_id=1";
  $result = mysql_query($query,$con);// or die("Bad query:$query");
  $numRows =mysql_numrows($result);
  $aa="";
  if($numRows>0){
     $row=mysql_fetch_assoc($result);
     $aa=$row["aa"];
  }
  return $aa;
}
function dna2rna($codon){
 $data=$codon;
 if($data[0]=="T")$data[0]="U";
 if($data[1]=="T")$data[1]="U";
 if($data[2]=="T")$data[2]="U";
 return $data;
}  
/*********************************************************
  this function displays all the   
  data (strains) by providers in a tabular format
***********************************************************/
function getStrainBySourceTable(){
  $strains =getStrainMainList(); global $SNPSOURCE;
   // display the columns header
  $source_name="Strain Name";
  $len=strlen($source_name);
  $verticalName=""; $straincountlist="<tr><td>Strains</td>";
  for($p=0;$p<=$len-1;++$p)$verticalName.="$source_name[$p]<br/>";
  $data="<table border='1' id='strains'><tr><th valign='top'>$verticalName</th>";
  $snp_sources=""; $i=0;
  while(list($source_id,$source_name)=each($SNPSOURCE)){
        $strain_count=getStrainCount($source_id);
        $len=strlen($source_name);
        $verticalName="";
        for($p=0;$p<=$len-1;++$p)$verticalName.="$source_name[$p]<br/>";
        $data.="<th class='strainlabel' valign='bottom' >$verticalName </th>";
        $straincountlist.="<td>$strain_count</td>";
        if($i==0)$snp_sources="$source_id";
        else $snp_sources .=":$source_id";
        ++$i;
  }
  $data.="</tr>";
  $data .="$straincountlist</tr>";
  $strainSources=split(":",$snp_sources); $i=0;
  // now for each strain, get strain_id,for e
  while(list($strain_id,$strain_name)=each($strains)){
     $data.="<tr><td class='sectionName' >$strain_name </td>";
     $i=0;
     while($i<count($strainSources)){
           $source_id=$strainSources[$i];
           $train_in_source=is_StrainInSource($strain_id,$source_id);
           $is_inlist="<td></td> ";
           if($train_in_source==1)$is_inlist="<td class='hasit'></td>";
            $data .="$is_inlist";
           ++$i;
      }
      $data.="</tr>";     
   }
   $data .="</table>";
   return $data;
}
/*********************************************************
  given a source_id and a strain_id, this function returns  
  true if this strain was provided by this source
***********************************************************/
function is_StrainInSource($strain_id,$snp_source){
    global $con;
    if(!$con)$con=getConnection();
    $query="select strain_id ";
    $query.=" from snp_strain_by_source";
    $query.=" where source_id=$snp_source and strain_id=$strain_id";
    $result = mysql_query($query,$con);// or die("Bad query $query:".mysql_error());
    $numRows = mysql_numrows($result);
    $train_in_source=0;
    if($numRows>0){
      $train_in_source=1;
    }
    return $train_in_source;
}

/**************************************************************
 given a source_id, returns the strain count
***************************************************************/
function getStrainCount($source_id){
   global $con;
   if(!$con)$con=getConnection();
    $query="select count(distinct strain_id) as strain_count ";
    $query.=" from snp_strain_by_source";
    $query.=" where source_id=$source_id ";
    $result = mysql_query($query,$con);
    $numRows = mysql_numrows($result);
    $train_count=0;
    if($numRows>0){
       $train_count=  mysql_result($result,$i,"strain_count");;
    }
    return $train_count;

}
/***************************************************************
   Given a SNP_id
   This function generates the list of strains with their IDs
****************************************************************/
function getStrainList($con,$snpid,$imputed,$selectedStrains){
        $strains="";  global $con;
        if(!$con)$con=getConnection();
        $strain_list=split(":",$selectedStrains);
       if(is_array($strain_list)&& $strain_list[0]){
        $query="select distinct s.strain_id,strain_name from snp_strain s,";
        if($imputed==0)$query .=" snp_genotype sg ";
        else
           $query.=" snp_imputed sg ";
        $query .=" where sg.snpid=$snpid and sg.strain_id=s.strain_id"; 
        $i=1;
       $query .=" and s.strain_id in ($strain_list[0]";
        while($i<count($strain_list)){
           $query .=",$strain_list[$i]";
           ++$i;
        }
        $query .=")";
        $query.="  order by s.strain_name";
        $result = mysql_query($query,$con);
        $numRows = mysql_numrows($result);
        $i=0;
        if($numRows>0)
        {
           while($i<$numRows){
                 $strainID= mysql_result($result,$i,"strain_id");
                 $strainName = mysql_result($result,$i,"strain_name");
                 $strains{$strainID}=$strainName;
                 $i++;
            }
        }
      }
     return $strains;
}


/*********************************************************
  given a codon, this function returns the frequency 
  of the the codon with max codon_fraction_1000 
***********************************************************/
function getCodonFrequency($codon,$con){
  $rnaCodon=dna2rna($codon);
  $query="select codon_fraction_1000 as codon_freq ";
  $query .=" from cgd_mouse_codonFrequency ";
  $query .=" where codon ='$rnaCodon' and organism_id=1";
  $result = mysql_query($query,$con);// or die("bad query:$query");
  $numRows =mysql_numrows($result);
  $codon_freq;
  if($numRows>0){
     $row=mysql_fetch_assoc($result);
     $codon_freq=$row["codon_freq"];
  }
  return $codon_freq;
}  
/****************************************************
* displayPages:
* Displays the number of pages of the results
* of a given and a link to each page  **************
*****************************************************/
function displayPages($totalSNPS,$filterlink,$rowcount,$page,$snp_loc,$fields){
  
   $numpages= ceil($totalSNPS/$rowcount);
   if(($page==0)||($page>$numpages))$page=1;
   $options=""; $type=2;
   if($fields{"v_type"}==1) $options.="&vt=1";
   if($fields{"g_name"}==1) $options.="&gn=1";
   if($fields{"f_class"}==1) $options.="&fc=1";
   if($fields{"cpg"}==1) $options.="&cp=1";
   if($fields{"m_allele"}==1) $options.="&ma=1";
   if($fields{"m_genotype"}==1) $options.="&mg=1";
   // compute next page number and previous page number
   $j=1;
  $data="<div style='width:90%;height:1em;margin:auto;padding:0.4em;font-size:0.75e;border-top:solid 1px #ddd;'>
         <center><b>Page(s): </b>";
   $prev=0;
   // re-ajust j. Found the index of the first page on this block
   // such that index=k-1
   $k=1;
   while($k<=$page){
       if($k%10==0){
          if(($page-$k)<10){$prev=$k-1;break;}
          elseif(($page-$k)==10){
              $k=$page;$prev=$page-1;break;}
       }
       ++$k;
   }
   if($page<10){$j=1;$prev=0;}
   else{$j=$k;}
   // display a link to the previous page
 if($prev>0){
    $pages .="<span class='pages'>";
    $pages.="<a name='page=$prev' 
             onclick=\"updatePageResult($type,this);\" >Previous<<</a></span>";
  }
  $data.="$pages";$pages="";
  while($j<=$numpages){
          if($j==$page){
             $pages.="<span id='pag' style='color:#804000;'>[$page]</span>";
          }
          else{
              $pages = "<span class='pages'>";
              $pages.="<a onclick=\"updatePageResult($type,this);\" name='page=$j' href='#'>";
             if((($j%10)==0)&&($j>$page)){
                  $pages .=">>Next</a></span>";
                  break;
              }
              else{ $pages .="[ $j</a>]</span>";}
           }
           $data.="$pages";
           $pages="";
          ++$j;
     }
     $data.="</font></center></div>";
     return $data;
 }
 
/*********************************************************
 getEnsGene:
 Given a transcrip_local_id, returns Ensembl gene_id
 **********************************************************/
function getEnsGene($gene_id){
    $gene_accession_id="";
    global $con;
    if(empty($con))$con=getConnection();
    $query="select accession_id from cgd_genes_desc  ";
    $query .=" where gene_id=$gene_id";
    $result = mysql_query($query,$con);
    $numRows = mysql_numrows($result);
    if($numRows){
       $row=mysql_fetch_assoc($result);
       $gene_accession_id=$row["accession_id"];
    }
    return $gene_accession_id;
}
/*********************************************************
 getMgiGene:
 Given a transcrip_local_id, returns MGI gene_id
 **********************************************************/
function getMgiGene($gene_id){
    global $con;
    if(empty($con))$con=getConnection();
    $query="select distinct mgi_geneid,marker_symbol from cgd_genes_ensembl_mgi ";
    $query .=" where  gene_id=$gene_id";
    $result = mysql_query($query,$con);
    $numRows =mysql_numrows($result);
    $mgi_gene="";
    if($numRows>0){
       $row=mysql_fetch_assoc($result);
       $mgi_gene[0]=$row["mgi_geneid"];
       $mgi_gene[1]=$row["marker_symbol"];
    }
    return $mgi_gene;
} 
/*********************************************************
 getEntrezGene:
 Given a gene_id, returns Entrez gene_id
 **********************************************************/
function getEntrezGene($gene_id){
    $query="select distinct entrez_geneid,marker_symbol from cgd_transcripts ts, cgd_genes_ensembl_entrez m";
    $query .=" where ts.gene_id=m.gene_id";
    if($gene_id>0)$query .=" and ts.gene_id=$gene_id ";
    $query .=" limit 1";
     global $con;
    if(empty($con))$con=getConnection();
    $result = mysql_query($query,$con);
    $numRows =mysql_numrows($result);
    $entrez_gene="";
    if($numRows>0){
       $row=mysql_fetch_assoc($result);
       $entrez_gene[0]=$row["entrez_geneid"];
       $entrez_gene[1]=$row["marker_symbol"];
    }
    return $entrez_gene;
} 

/*********************************************************
 getTranscriptInfo:
 Given a transcript_id, returns Ensembl transc accession_id
 transcript_start, transcript_end, transcript_strand,
 **********************************************************/
function  getTranscriptInfo($transcript_id){
   $query =" select accession_id,transcript_start, transcript_end, transcript_strand ";
   $query .=" from cgd_transcripts ts, cgd_transcripts_desc td ";
   $query .=" where ts.transcript_local_id= td.transcript_local_id ";
   if($transcript_id>0) $query .=" and ts.transcript_local_id=$transcript_id";
   $query .=" limit 1";
   global $con;
   if(empty($con))$con=getConnection();
   $result = mysql_query($query,$con); 
   $accession_id=""; $transcdata=""; $strand="+";
   $numRows = mysql_numrows($result);
   if($numRows>0){
      $accession_id=mysql_result($result,0,"accession_id");
      $tstart= mysql_result($result,0,"transcript_start");
      $tend= mysql_result($result,0,"transcript_end");
      $tstrand= mysql_result($result,0,"transcript_strand");
      $transcdata{"id"}=$accession_id; $transcdata{"tstart"}=$tstart; $transcdata{"tend"}=$tend; 
      $transcdata{"strand"}=$tstrand;
    }
    return $transcdata;
   
}
/*************************************************************
  returns the number of exons for a given transcript
**************************************************************/
function getExonCount($transcript_local_id){
       // get the exon rank
       $query="select count(*) as exon_count from cgd_transcript_exons  where transcript_local_id=$transcript_local_id ";         
       global $con;
       if(empty($con))$con=getConnection();
       $result = mysql_query($query,$con); 
       $Rows = mysql_numrows($result);
       $exon_count=0;
       if($Rows>0){$exon_count=mysql_result($result,0,"exon_count");}
       return $exon_count;
}
/*************************************************************
  returns the exon rank of a given exon for a given transcript
**************************************************************/
function getExonRank($transcript_local_id,$bp_position){
       // get the exon rank
       $query="select exon_rank from cgd_transcript_exons  where transcript_local_id=$transcript_local_id ";
       $query .=" and exon_start<=$bp_position and exon_end>=$bp_position";                  
       global $con;
       if(empty($con))$con=getConnection();
       $result = mysql_query($query,$con); 
       $Rows = mysql_numrows($result);
       $exon_rank=0;
       if($Rows>0){$exon_rank=mysql_result($result,0,"exon_rank");}
       return $exon_rank;
}
/*************************************************************
  returns the exon rank of a given exon for a given transcript
**************************************************************/
function getIntronRank($transcript_local_id,$bp_position){
       // get the exon rank
       $query="select et.exon_rank as ex1_rank,et2.exon_rank as ex2_rank "; //, et2.exon_rank as ex2_rank ";
       $query.="from cgd_transcript_exons et,cgd_transcript_exons et2 ";
       $query.=" where et.transcript_local_id=$transcript_local_id ";
       $query .=" and et.transcript_local_id=et2.transcript_local_id ";
       $query .=" and et.exon_end< $bp_position and  et.exon_id!= et2.exon_id ";
       $query .=" and ((et.exon_rank+1=et2.exon_rank)or (et.exon_rank=et2.exon_rank+1))";
       $query.="  and et2.exon_start >$bp_position ";                
       global $con;
       if(empty($con))$con=getConnection();
       $result = mysql_query($query,$con); 
       $Rows = mysql_numrows($result);
       $intron_rank=0;
       if($Rows>0){
          $intron_rank=mysql_result($result,0,"ex1_rank");
          $ex2_rank=mysql_result($result,0,"ex2_rank");
          if($intron_rank>$ex2_rank)
             $intron_rank=$ex2_rank;
        }
       return $intron_rank;
}
function getProteinData($transcript_local_id){
 $query="select cds_start,cds_end,protein_id from cgd_ensembl_protein  p";
 $query .=" where p.transcript_local_id=$transcript_local_id";
  global $con;
  if(empty($con))$con=getConnection();
  $result = mysql_query($query,$con); //or die ("No result for your query:$query");
   $Rows = mysql_numrows($result);
   $translation="";
   if($Rows>0){
      $cds_start=mysql_result($result,0,"cds_start");
      $cds_end= mysql_result($result,0,"cds_end");
      $protein= mysql_result($result,0,"protein_id");
      $translation="$cds_start:$cds_end:$protein";
   }
  return $translation;

}
/*********************************************************************
 returns the upstream and downstream gene information for a given
 intergenic SNP. Gene should be within +- 5kb from the SNP position
**********************************************************************/
function getIntergenic($chromosome_id,$bp_position){
    global  $MGI_GENEDETAIL;
     $query="select distinct gd.accession_id as geneid,e.gene_start,e.gene_end,gd.gene_type,ts.transcript_strand,ts.gene_id";
     $query .="  from cgd_genes e,(select chromosome_id,max(gene_end) as gene_end2 from cgd_genes where chromosome_id=$chromosome_id ";
     $query .=" and gene_end<=$bp_position  group by chromosome_id) as t,cgd_genes_desc gd,cgd_transcripts ts";
     $query .="  where e.chromosome_id =t.chromosome_id";
     $query .=" and e.gene_end =t.gene_end2";
     $query .=" and e.gene_id =ts.gene_id";
     $query .=" and e.gene_id =gd.gene_id";
     global $con;
     if(empty($con))$con=getConnection();
     $result = mysql_query($query,$con);
     $numRows = mysql_numrows($result);
     $i=0;
     $data="<table border='1'><caption>Nearest Genes</caption>
             <tr class='sectionOptions'><th class='innertranstable'>Ensembl Gene_id</th><th class='innertranstable'>Gene Start</th>";
     $data.="<th class='innertranstable'>Gene End</th><th class='innertranstable'>Gene Type</th>";
     $data.=" <th class='innertranstable'>Orientation</th><th class='innertranstable'>SNP_Position</th>
              <th class='innertranstable'>Distance from SNP</tr>";
     if($numRows>0){ //display all genes
        while($i<$numRows){
                   $geneid=mysql_result($result,$i,"geneid");
                   $gene_start=mysql_result($result,$i,"gene_start");
                   $gene_end=mysql_result($result,$i,"gene_end");
                   $gene_type=mysql_result($result,$i,"gene_type");
                   $gene_strand=mysql_result($result,$i,"transcript_strand");
                   $distance=$bp_position-$gene_end;
                   $gene_id=mysql_result($result,$i,"gene_id");
                    if($gene_strand=="1")$gene_strand="+";
                    // for each gene, get the mgi symbol
                   $query="select distinct mgi_geneid,marker_symbol from cgd_transcripts ts, cgd_genes_ensembl_mgi m,";
                   $query .=" (select gene_id,max(transcript_end)as max_end from cgd_transcripts group by gene_id)as t ";
                   $query .=" where  ts.gene_id=$gene_id and ts.gene_id=t.gene_id and ts.transcript_end=t.max_end ";
                   $query .=" and ts.gene_id=m.gene_id";
                   $mgigene="";
                   $mgigenename="";
                   $result2 = mysql_query($query); //or die ("No result for your query: $query");
                   $numRows2 = mysql_numrows($result2);
                   if($numRows2>0){
                       $mgigeneid=mysql_result($result2,0,"mgi_geneid");
                       $mgigenename=mysql_result($result2,0,"marker_symbol");
                       $mgigene="<br> MGI Symbol: <a href='$MGI_GENEDETAIL/accession_report.cgi?id=$mgigeneid'> $mgigenename</a>";
                   }
                   $data.="<tr class='sectionOptions' style='background:#fff;'><td>$geneid $mgigene";
                   $data.=" <!-- <a href=\"http://cgd.jax.org/cgdsnpdb/search/?geneViw=$geneid\"> [GeneSnpsView]</a>--></td>";
                   $data.= "<td > $gene_start</td>";
                   $data.="<td > $gene_end</td><td >$gene_type</td>";
                   $data.=" <td > $gene_strand</td><td> $bp_position</td>";
                   $data.=" <td>$distance bp (Left)</td></tr>";
                  $i++;
              }
         }
          // now get the downstream
           $query ="select distinct gd.accession_id as geneid,e.gene_start,e.gene_end,gd.gene_type,ts.transcript_strand,ts.gene_id";
           $query .="  from cgd_genes e,(select chromosome_id,min(gene_start) as gene_start 
                       from cgd_genes where chromosome_id=$chromosome_id ";
           $query .=" and gene_start>=$bp_position  group by chromosome_id) as t,cgd_genes_desc gd,cgd_transcripts ts";
           $query .="  where e.chromosome_id =t.chromosome_id";
           $query .=" and e.gene_start =t.gene_start";
           $query .=" and e.gene_id =ts.gene_id";
           $query .=" and e.gene_id =gd.gene_id";
           $result = mysql_query($query); //or die ("No result for your query");
           $numRows = mysql_numrows($result);
            $i=0;
           if($numRows>0){ //display all genes
              while($i<$numRows){
                   $geneid=mysql_result($result,$i,"geneid");
                   $gene_start=mysql_result($result,$i,"gene_start");
                   $gene_end=mysql_result($result,$i,"gene_end");
                   $gene_type=mysql_result($result,$i,"gene_type");
                   $gene_strand=mysql_result($result,$i,"transcript_strand");
                   $gene_id=mysql_result($result,$i,"gene_id");
                   $distance=$gene_start-$bp_position;
                   if($gene_strand=="1")$gene_strand="+";
                   // for each gene, get the mgi symbol
                   $query="select distinct mgi_geneid,marker_symbol from cgd_transcripts ts, cgd_genes_ensembl_mgi m,";
                   $query .=" (select gene_id,min(transcript_start)as min_start from cgd_transcripts group by gene_id)as t ";
                   $query .=" where  ts.gene_id=$gene_id and ts.gene_id=t.gene_id and ts.transcript_start=t.min_start ";
                    $query .=" and ts.gene_id=m.gene_id";
                   $mgigene="";
                   $mgigenename="";
                   $result2 = mysql_query($query); //or die ("No result for your query");
                   $numRows = mysql_numrows($result2);
                   if($numRows>0){
                       $mgigeneid=mysql_result($result2,0,"mgi_geneid");
                       $mgigenename=mysql_result($result2,0,"marker_symbol");
                        $mgigene="<br> MGI Symbol: <a href='$MGI_GENEDETAIL/accession_report.cgi?id=$mgigeneid'> $mgigenename</a>";
                   }
                   $data.="<tr class='sectionOptions' style='background:#eee;' ><td >$geneid $mgigene";
                   $data.=" <!-- <a href=\"http://cgd.jax.org/cgdsnpdb/search/?geneViw=$geneid\"> [GeneSnpsView]</a>--></td>
                           <td > ";
                   $data.="$gene_start</td><td> $gene_end</td><td >$gene_type</td>";
                   $data.=" <td> $gene_strand</td><td> $bp_position</td>";
                   $data.= "<td >$distance bp (Right)</td></tr>";
                  $i++;
               }
           }
          $data.="</table>";
         return $data;
}
/*********************************************************
 getSNPAccessionSource:
 Given a SNP accession id, returns the source id list
 **********************************************************/
function getSNPAccessionSource($snpAccession){
   // get the cgd snpid of this accession_id
   $snpid=getSNPID($snpAccession);
   $query =" select distinct source_id from snp_by_source ";
   $query .=" where snpid=$snpid";
   global $con;
   if(empty($con))$con=getConnection();
   $result = mysql_query($query,$con); 
   $source_id=0;$i=0;
   $numRows = mysql_numrows($result);
   if($numRows>0){
      while($i<$numRows){
            $id=mysql_result($result,$i,"source_id");
            if($i==0)$source_id="$id";
            else{$source_id.=":$id";}
            ++$i;
       }
   }
   return $source_id;
   
}
/*********************************************************
 getSNPID:
 Given a SNP accession id, returns the cgd snpid
 **********************************************************/
function getSNPID($snpAccession){
   $query =" select snpid from snp_accession ";
   $query .=" where accession_id='$snpAccession'";
   global $con;
   if(empty($con))$con=getConnection();
   $result = mysql_query($query,$con); 
   $snpid=0;$i=0;
   $numRows = mysql_numrows($result);
   if($numRows>0){
      $snpid=mysql_result($result,0,"snpid");
    }
    return $snpid;
   
}
/******************************************************************
 returns true if a given snp was provided by the usere 
 selected snp sources list
*******************************************************************/
function is_snpInSelectedSource($selectedSources,$thisSNP_sources){
   $sources= split(":",$selectedSources); $found=0;
   foreach($sources as $source_id){
          if($thisSNP_sources{"$source_id"}==$source_id){
                ++$found; break;
          }
    }
    return $found;
}
/*************************************************************************** 
   getSessionQuerySNPsCount: given a genomic region,source(s),and annotation(s)
    this function returns a count of SNPs Ids within the specified region
*****************************************************************************/
function getSessionQuerySNPsCount($chr,$startpos,$endpos,$selectedSources,$snp_location){
   $snplist=$_SESSION['genome']["$chr"]["$startpos"]["$endpos"]["$snp_location"]["snpList"];
   //print_r($snplist);
   $count=0;
   foreach($snplist as $snpid){
          //get this snp sources
         $thisSNP_sources= $_SESSION['snp']["$snpid"]['source'];
         if(is_snpInSelectedSource($selectedSources,$thisSNP_sources)){
               if($snp_location==12)++$count;
               else{ //check both the source and snp_location match
                   $thisSNP_loc= $_SESSION['snp']["$snpid"]['loc'];
                   $snploc=array();
                   if(($snp_location!=10)&&($snp_location!=2))$snploc[0]=$snp_location;
                   else if($snp_location==10)$snploc=array(1,2,3,4,5,6,7,8);
                   else if($snp_location==2)$snploc=array(2,5,6,7,8);
                   foreach($snploc as $loc_id){
                      if(($loc_id>0)&&($thisSNP_loc{"$loc_id"}==$loc_id)){
                         ++$count; break;
                      }
                   }
               }
         }
   }return $count;
}
/*************************************************************************** 
   getQuerySNPsCount: given a genomic region,source(s),and annotation(s)
   this function returns a count of SNPs Ids within the specified region
*****************************************************************************/
function getQuerySNPsCount($con,$chr,$startpos,$endpos,$selectedSources,$snp_location,$gene_id,$transcript_local_id,$sg){
   global $con; $count=0;if(!$con)$con=getConnection();
    // get the user selected sources
   $sources= split(":",$selectedSources); 
   $query="select count(distinct sb.snpid) as snpcount ";
   $query .=" from snp_by_source sb , snp_position m ";
   if(($gene_id>0)|| ($transcript_local_id>0)){ //query using gene annotations
         $query .=" ,snp_transcript st where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos";
         $term=" st.transcript_local_id=$transcript_local_id ";
         if($gene_id>0)$query.=" and st.gene_id=$gene_id "; if($transcript_local_id>0)$query .=" and $term ";
        //check the location 
        if(($snp_location!=-1)&& ($snp_location!=12)){
            if($snp_location==10){$query .=" and st._loc_func_key in(1,2,3,4,5,6,7,8)";}
            else if($snp_location==2){$query .=" and st._loc_func_key in(2,5,6,7,8)";}
            else if($snp_location==34){$query .=" and st._loc_func_key in(3,4)";}
            else $query.=" and st._loc_func_key =$snp_location ";
        }$query .="  and st.snpid=sb.snpid";
    }
   else{
       if(($snp_location==-1)||($snp_location==12)){
          $query.=" where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
          if($snp_location==-1)$query .=" and m.is_intergenic=1 ";
        }
       else{
          if(($snp_location==-1)){
              $query.=", snp_main ms where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
              $query .=" and m.snpid=ms.snpid and ms.is_intergenic=1 ";
          }
          else if(($snp_location==12)){ // all snps
             $query.=" where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos";
          }
         else{
             $query.=",snp_transcript st where m.chromosome_id=$chr and m.bp_position between $startpos and $endpos ";
             $query .=" and m.snpid=st.snpid ";
             if($snp_location==10){$query .=" and st._loc_func_key in(1,2,3,4,5,6,7,8)";}
             else if($snp_location==2){$query .=" and st._loc_func_key in(2,5,6,7,8)";}
             else if($snp_location==34){$query .=" and st._loc_func_key in(3,4)";}
             else $query.=" and st._loc_func_key =$snp_location ";
          }
       } 
    }
    $query .=" and sb.snpid=m.snpid  and sb.source_id in(";
    $i=0;
    foreach($sources as $source_id){
           if($i==0)$query .="$source_id";
           else{
              if($source_id>0)$query .=",$source_id";
            }
          ++$i;
    }
    $query .=") ";
    $result = mysql_query($query,$con);
    $numRows = mysql_numrows($result); $i=1;
    $row=mysql_fetch_assoc($result);
    $count=$row["snpcount"];
    return $count;
}
/*********************************************************************************
   getMinorAlleleFreq returns an array where index0=strain_count,
   index1=Ns_count,index2=minor_allele_count
   index3= minor_allele for a given snp and source(s)
**********************************************************************************/
function getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$imputed){ 
        global $con;if(!$con)$con=getConnection();$accession="";$genotypeSum="";
        $sources=split(":",$selectedSources);
        if($imputed==0){
            // $query ="select distinct strain_id,genotype_allele from snp_genotype where snpid=$snpid and source_id in (";
            $query ="select source_id,genotype_allele, count(distinct strain_id) as strain_count ";
            $query.=" from snp_genotype where snpid=$snpid and source_id in (";
            $i=0;
            foreach($sources as $source_id){
               if($source_id>0){
                  if($i==0)$query.="$source_id";
                  else{$query.=",$source_id";}
                  ++$i;
                }
            }
           $query .=") group by source_id,genotype_allele" ; //order by strain_id";
        }
        else{  
           $query ="select genotype_allele, count(distinct strain_id) as strain_count  from snp_imputed where snpid=$snpid ";
           $query.=" group by genotype_allele"; //order by strain_id";
        }
        $result = mysql_query($query);// or die ("Bad query:$query");
        $numRows = mysql_numrows($result);
        $i=0;$strain_count=0;$minor_allele="";
        $minor_allele_count=0;$Ns_count=0; 
        $b6_allele_count=0;$next_strain=0;
        if($numRows>0){
           if($imputed==0){
              while($i<$numRows){
                   $source_id=mysql_result($result,$i,"source_id");
                   $strain_count=mysql_result($result,$i,"strain_count");
                   $genotype_allele= mysql_result($result,$i,"genotype_allele");
                   $genotypeSum{$source_id}{$genotype_allele}=$strain_count;
                   //print "--$snpid-$source_id-$genotype_allele--";
                   ++$i;
               }
           }
           else{
              $source_id=21;
              while($i<$numRows){
                   $strain_count=mysql_result($result,$i,"strain_count");
                   $genotype_allele= mysql_result($result,$i,"genotype_allele");
                   $genotypeSum{$source_id}{$genotype_allele}=$strain_count;
                   //print "--$snpid-$source_id-$genotype_allele--$strain_count -=";
                   ++$i;
               }
           }
        }
       return $genotypeSum;
}
/***********************************************************************
  for a given source_id,and snp, this fuction returns a map
  of b6_allele count,snp_allele_count,N_count,H_count, and ambigo=uous count
************************************************************************/
function parseMinorAllele($genoSumarry,$impgenoSum,$snp_allele,$black6_allele){
   $genomap="";
   if(is_array($genoSumarry)){
      while(list($source_id,$genoallele)=each($genoSumarry)){                
              $totalsourceStrain=0;$genomap{$source_id}{"ambiguous"}=0;$genomap{"H_allelecount"}{$source_id}=0;
              $genomap{"N_allelecount"}{$source_id}=0;
              while(list($geno_allele,$straincount)=each($genoallele)){
                      //$minorallele.="$source_id:$geno_allele:$straincount -";  
                      $totalsourceStrain+=$straincount; 
                       if(strtoupper($geno_allele)==strtoupper($black6_allele))
                           $genomap{"b6_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)==strtoupper($snp_allele))
                           $genomap{"snp_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)=="N")
                           $genomap{"N_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)=="H")
                           $genomap{"H_allelecount"}{$source_id}=$straincount;
                        else
                          $genomap{"ambiguous"}{$source_id}=$straincount;
              }
             // $minorallele.="$source_id:$genoallele:$straincount -"; 
             $genomap{"straincount"}{$source_id}=$totalsourceStrain;
        }
     }
     if(!empty($impgenoSum)){
               $source_id=16;
               $totalsourceStrain=0;$genomap{"ambiguous"}{$source_id}=0;$genomap{"H_allelecount"}{$source_id}=0;
               $genomap{"N_allelecount"}{$source_id}=0;
               $genoallele=$impgenoSum{$source_id};
               while(list($geno_allele,$straincount)=each($genoallele)){
                        //$minorallele.="$source_id:$geno_allele:$straincount -"; 
                        $totalsourceStrain+=$straincount; 
                        if(strtoupper($geno_allele)==strtoupper($black6_allele))
                           $genomap{"b6_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)==strtoupper($snp_allele))
                           $genomap{"snp_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)=="N")
                           $genomap{"N_allelecount"}{$source_id}=$straincount;
                        else if(strtoupper($geno_allele)=="H")
                           $genomap{"H_allelecount"}{$source_id}=$straincount;
                        else
                          $genomap{"ambiguous"}{$source_id}=$straincount;
               }
              $genomap{"straincount"}{$source_id}=$totalsourceStrain;
     }
     return $genomap;
            
}
/*********************************************************************
 returns the the minor allele info including the allele, frequency,
 strain count
**********************************************************************/
function computeMAF($genomap,$black6_allele,$snp_allele,$sourcelist){
  $Ns;$minMaf=0;$maxMaf=0;$minal=""; $minmap="";$nmap="";
  $t=0;$minalt=""; $minmapt="";$nmapt=""; $totalStrain="";
 //print_r($genomap);
  foreach($sourcelist as $source_id){
          $Ns=$genomap{"N_allelecount"}{$source_id};$strcount=$genomap{"straincount"}{$source_id};
          $amb=$genomap{"ambiguous"}{$source_id};$Hs=$genomap{"H_allelecount"}{$source_id};$genocount=0;
          $ballelecount=$genomap{"b6_allelecount"}{$source_id}; 
          $snpallelecount=$genomap{"snp_allelecount"}{$source_id};
                 //get the b6 count
          if($genomap{"b6_allelecount"}{$source_id}<=$genomap{"snp_allelecount"}{$source_id}){
             $minalt=$black6_allele;                      
             $genocount=$genomap{"b6_allelecount"}{$source_id}; 
          }
          else{
             $minalt=$snp_allele;                      
             $genocount=$genomap{"snp_allelecount"}{$source_id}; 
          }
          if(($strcount-$Ns -$Hs)>0){
             $minmaft= ceil(($genocount/($strcount-$Ns -$Hs))*100);
             $nmapt =ceil(($Ns/$strcount)*100);
          }
          if(($ballelecount>0)&&($snpallelecount>0)){
              if($t==0){
                 $minal=$minalt;
                 $minmaf ="$minmaft ";
                 $nmap =$nmapt; ++$t;
                 $totalStrain="$strcount";
              }
              else{
                 $minal  .=":$minalt";
                 $minmaf .=":$minmaft ";
                 $nmap   .= ":$nmapt";
                 $totalStrain.=":$strcount";
              }
          }
   }
   $data{"minal"}=$minal;
   $data{"minmaf"}=$minmaf;
   $data{"nmap"}=$nmap;
   $data{"totalStrain"}=$totalStrain;
   return $data;
}

/*********************************************************
getSNPConflictList:
Given a snpid, returns all the issues associsted with it
**********************************************************/
function  getSNPConflictList($con,$snpid){
   $query =" select distinct error_flag from snp_error_log ";
   $query .=" where snpid=$snpid"; 
   global $con;
   if(!$con)$con=getConnection();
   $result = mysql_query($query,$con); //or die ("Bad query:$query"); 
   $numRows = mysql_numrows($result);
   $conflictlist="";$i=0;
    while($i<$numRows){
      $error_type=mysql_result($result,$i,"error_flag");
      if($error_type==6)$conflictlist.="<li>Conflicting SNP allele between sources</li>";
      if($error_type==7)$conflictlist.="<li>Conflicting strains genotype between SNP sources</li>";
      if($error_type==9)$conflictlist.="<li>Ambiguous SNP due to multiple genome alignments</li>";
      if($error_type==10)$conflictlist.="<li>SNP located in ambiguous cluster window</li>";
       $i++;
    }
    return $conflictlist;
   
}
/*********************************************************

Given a snpid, returns all the sources associa`ted with it
**********************************************************/
function getSNPsources($con, $snpid){
       global $con;
       if(!$con)$con=getConnection();
      $query="select s.source_id,s.source_name,s.website ";
      $query .=" from snp_source s,snp_by_source sb";
      $query .="  where sb.snpid =$snpid";
      $query .="  and sb.source_id =s.source_id";
      $result = mysql_query($query,$con); //or die ("No result for your query");
      $numRows = mysql_numrows($result);
      $i=0;$source_list="";$source_name="";$source_link="";$imputed=0;$sources="";$source_ids="";
      if($numRows>0){ // snp has source
         while($i<$numRows){
                $source_name=mysql_result($result,$i,"source_name");
                $source_id= mysql_result($result,$i,"source_id");
                $source_link= mysql_result($result,$i,"website");
                
                if($i==0){
                   $source_ids="$source_id";
                   $source_list="<a href='$source_link' target='_blank'> $source_name </a>";
                }
                else{
                   $source_ids.=":$source_id";
                   $source_list.="&nbsp;||&nbsp;<a href='$source_link' target='_blank'> $source_name </a>";
                }
                if(($source_id==16)||($source_id==21))$imputed=1;
                $i++;
         }
         $sources{"list"}=$source_list;
         $sources{"type"}=$imputed;
         $sources{"ids"}=$source_ids;
      }
      return $sources;
}
/*************************************************************
 Given the snpid and transcriptid,
 this function returns the polarity and the hydropathy
 index information of the aminoacid coded by this SNP
***************************************************************/
function getAaPolarity($snpid,$transcript_local_id,$type){
    $data;
    if($snpid>0 && $transcript_local_id>0){
       global $con;
       if(!$con)$con=getConnection();
       $query="select ac.acidity_polarity as b6_aa_polarity,ac.hydropathy_index as b6_hydropathy";
       $query .=" from snp_aminoacid a,cgd_aminoacid ac";
       $query .="  where snpid=$snpid  and a.transcript_local_id=$transcript_local_id ";
       if($type==1){
          $query .=" and a.ref_aa=ac.three_letter_abrev ";
       }
       else{
          $query .=" and a.snp_aa=ac.three_letter_abrev ";
       }
       $result = mysql_query($query,$con);//or die ("Invalid :");
       $Rows = mysql_numrows($result);
       if($Rows>0){
          $b6_aa_polarity=mysql_result($result,0,"b6_aa_polarity");
          $b6_hydropathy=mysql_result($result,0,"b6_hydropathy");
          $data{"polarity"}=$b6_aa_polarity;
          $data{"hydropathyindex"}=$b6_hydropathy;
       }
    }
   return $data;
}
/***************************************************************
  Given the aminoacid three letter symbol,
  this function returns the equivalent one letter abreviation
****************************************************************/
function getAaOneletter($aa){
  $oneletter="";
  if(!empty($aa)){
      global $con;
      if(!$con)$con=getConnection();
      $query="select one_letter_abrev from cgd_aminoacid where three_letter_abrev='$aa'";
      $result = mysql_query($query,$con); //or die ("Invalid query");
      $Rows = mysql_numrows($result);
      if($Rows>0){ 
         $b6leter=mysql_result($result,0,"one_letter_abrev");
         $oneletter{"char"}=$b6leter;
         $oneletter{"str"}="($b6leter)";
       }
   }
  return $oneletter;
}
/*******************************************************************
  Given two aminoacid one letter symbols,
  this function returns the blosum62 score of the substitution
********************************************************************/
function getBlossum62Score($b6leter,$snpleter){
         $blossumScore=-99;
   if(!empty($snpleter)&&!empty($b6leter)){
         global $con;
         if(!$con)$con=getConnection();
         $query="select score";
         $query .=" from cgd_blosum62NCBI ";
         $query .="  where (amino_1='$b6leter' and amino_2='$snpleter') ";
         //$query .=" or (amino_1='$snpleter' and amino_2='$b6leter')";
         $result = mysql_query($query,$con); //or die ("Invalid query: $query");
         $Rows = mysql_numrows($result);
         if($Rows>0){
            $blossumScore=mysql_result($result,0,"score");
         }
    }
  return $blossumScore;
}
/*************************************************************************
   returns a table with a given substitution information
**************************************************************************/
function getTable($exon_rank,$frame,$snposincds,$snposinprotein,$blossumScore,$cai,$fb6_codon,$balck6_aa,$b6oneletter,$b6_codon,
                  $b6_aa_polarity,$b6_hydropathy,$b6hydropathy,$balck6_aa,$SNP_HOME2,
                  $fsnp_codon,$snp_aa,$oneletter,$snp_codon,$snp_aa_polarity,$snp_hydropathy,$snphydropathy,
                  $snp_aa){
     $locationtype=" <table border='1' style='background:#eee; font-size:0.85em;width:140;padding:0;margin:0;text-align:center;'>
                    <tr><td>Exon Rank</td></tr>
                    <tr><td style='background:#fff;'><b>$exon_rank</b></td><tr><td> SNP Codon Position</td></tr>
                    <tr><td style='background:#fff;'><b>$frame</b></td></tr>
                    <tr><td> SNP Position In CDS</td></tr>
                    <tr><td style='background:#fff;'><b>$snposincds</b></td></tr>
                    <tr><td> Codon Position In Protein</td></tr>
                    <tr><td style='background:#fff;'><b>$snposinprotein</b></td></tr>

                    <tr><td>BLOSUM62 replacement</td></tr>
                    <tr><td style='background:#fff;'><b>$blossumScore</b></td></tr>
                    <tr><td>&#916;CAI</td></tr>
                   <tr><td style='background:#fff;'><b>$cai</b></td></tr>
                  </table></td>
               <td class='sectionOptions'>
                   <table border='1'><caption>C57BL/6J Amino Acid</caption>
                       <tr><td  class='sectionName'>Codon Usage</td><td class='sectionOptions'>$fb6_codon</td></tr>
                       <tr><td  class='sectionName'>Aminoacid Symbol</td>
                           <td class='sectionOptions'>$balck6_aa$b6oneletter</td></tr>
                       <tr><td  class='sectionName'>Codon</td><td class='sectionOptions'>$b6_codon</td></tr>
                      <tr><td  class='sectionName'>Chemical Characteristics</td><td class='sectionOptions'>$b6_aa_polarity</td></tr>
                      <tr><td  class='sectionName'>Hydropathy Index</td>
                      <td class='sectionOptions'> $b6_hydropathy $b6hydropathy</td></tr>";
                        $b6_img=$balck6_aa;$b6_img .=".gif"; 
         $locationtype.="<tr><td class='sectionName'>Structures</td>
                              <td><img src='$SNP_HOME2/images/aminoacid/$b6_img' alt='$balck6_aa structure'/></td></tr>";
         $locationtype.="</table></td>";
         $locationtype.="<td class='sectionOptions'> 
                       <table border='1'><caption>SNP Variant Amino Acid</caption>
                            <tr><td  class='sectionName'>Codon Usage</td><td class='sectionOptions'>$fsnp_codon</td></tr>
                            <tr><td  class='sectionName'>Aminoacid Symbol</td><td class='sectionOptions'>$snp_aa$oneletter</td></tr>";
         $locationtype.="       <tr><td  class='sectionName'>Codon</td><td class='sectionOptions'>$snp_codon</td></tr>";
         $locationtype.="  <tr><td  class='sectionName'>Chemical Characteristics</td>
                                 <td class='sectionOptions'>$snp_aa_polarity</td></tr>";
         $locationtype.="  <tr><td  class='sectionName'>Hydropathy Index</td>
                               <td class='sectionOptions'> $snp_hydropathy $snphydropathy</td></tr>";
         $snp_img=$snp_aa;$snp_img .=".gif";
         $locationtype.="  <tr><td  class='sectionName'>Structures</td><td>
                                 <img src='$SNP_HOME2/images/aminoacid/$snp_img' alt='$snp_aa Structure'/></td></tr>";
          $locationtype.="  </table>";
     return $locationtype; 
}
/******************************************************************
  retuns all functional implications associated with a given snp
  on a specific transcript
******************************************************************/
function getSNPimplication($link,$snpid,$transcript_local_id,$bp_position){
   global  $SNP_HOME2; global $SNPCODON;
     $locationtype ="<br/><table bgcolor='#eeeee' border=1>";
     $query="select _frame_key,ref_aa,ref_codon,snp_aa,snp_codon,PosInCDS,PosInProtein from snp_aminoacid where snpid=$snpid 
             and transcript_local_id=$transcript_local_id ";
     global $con;
     if(!$con)$con=getConnection();
     $result = mysql_query($query,$con); //or die ("No result");
     $Rows = mysql_numrows($result);
     if($Rows>0){ // process cds snp
         $frame=mysql_result($result,0,"_frame_key");
          $balck6_aa=mysql_result($result,0,"ref_aa");
          $snp_aa=mysql_result($result,0,"snp_aa");
          $b6_codon=mysql_result($result,0,"ref_codon");
          // $snposincds,$snposinprotein
          $snposincds=mysql_result($result,0,"PosInCDS");
          $snposinprotein=mysql_result($result,0,"PosInProtein");
          $snp_codon=mysql_result($result,0,"snp_codon");
          $SNPCODON[0]=4;$SNPCODON[1]=$balck6_aa;$SNPCODON[2]=$b6_codon;$SNPCODON[3]=$snp_aa;$SNPCODON[4]=$snp_codon;
          $fb6_codon=getCodonFrequency($b6_codon,$con);
          $fsnp_codon=getCodonFrequency($snp_codon,$con);
          $fb6_max=getFmax($b6_codon,$con);
          $fsnpmax=getFmax($snp_codon,$con);
          $cai=round(abs(($fb6_max/$fb6_codon)*($fsnp_codon/$fsnpmax)),2);
          if($location==1){$locationtype .= "<tr bgcolor='#B0E0E6'>";}
          else{$locationtype .="<tr bgcolor='#F0E68C'>";}
          $snp_hydropathy=-99;$b6_hydropathy=-99;
          if($frame==1){
             $b6_codon= "<b><font color='red'>$b6_codon[0]</font>$b6_codon[1]$b6_codon[2]</b>";
             $snp_codon= "<b><font color='red'>$snp_codon[0]</font>$snp_codon[1]$snp_codon[2]</b>";
          }
          if($frame==2){
             $b6_codon= "<b>$b6_codon[0]<font color='red'>$b6_codon[1]</font>$b6_codon[2]</b>";
             $snp_codon= "<b>$snp_codon[0]<font color='red'>$snp_codon[1]</font>$snp_codon[2]</b>";
          }
          if($frame==3){
             $b6_codon= "<b>$b6_codon[0]$b6_codon[1]<font color='red'>$b6_codon[2]</font></b>";
             $snp_codon= "<b>$snp_codon[0]$snp_codon[1]<font color='red'>$snp_codon[2]</font></b>";
          }
          $snp_aa_polarity="Not Applicable";$b6_aa_polarity="Not applicable";
          // get chimical characteristics of the aminoacids of b6_aa
          if($balck6_aa!="stop"){
             $type=1;$polarity=getAaPolarity($snpid,$transcript_local_id,$type);
             $b6_aa_polarity=$polarity{"polarity"};$b6_hydropathy=$polarity{"hydropathyindex"};
          }
          if($snp_aa!="stop"){
             $type=0;$polarity=getAaPolarity($snpid,$transcript_local_id,$type);
             $snp_aa_polarity=$polarity{"polarity"};$snp_hydropathy=$polarity{"hydropathyindex"};
          }
          $b6oneletter="";$b6hydropathy="(Hydrophobic)";$snphydropathy="(Hydrophobic)";
          if($b6_hydropathy ==-99){
             $b6hydropathy="";$b6_hydropathy="Not Applicable";
          }
          else if($b6_hydropathy <0){$b6hydropathy="(Hydrophilic)";}
          if( $snp_hydropathy ==-99){
              $snphydropathy="";$snp_hydropathy="Not Applicable";
          }
          else if( $snp_hydropathy <0){$snphydropathy="(Hydrophilic)";}
          $b6leter;$snpleter;
          if($b6stop!=1){
               // get the one letter abreviation of this aminoacd
               $data=getAaOneletter($balck6_aa);$b6oneletter=$data{"str"};$b6leter=$data{"char"};
          }
          $oneletter="";
          if($snpstop!=1){  
               // get the one letter abreviation of this aminoacd
              $data=getAaOneletter($snp_aa);$oneletter=$data{"str"};$snpleter=$data{"char"};
          }
          $blossumScore=-99;
           // get the blossum62 matrix score of this substitution
           $snpleter=trim($snpleter); $b6leter=trim($b6leter);
          if(!empty($b6leter)&& !empty($b6leter)){
              $blossumScore=getBlossum62Score($b6leter,$snpleter);
           }
          $conservation="Not Applicable";
          if($blossumScore!=-99)$conservation="$blossumScore";
          $exon_rank=getExonRank($transcript_local_id,$bp_position);
          $locationtype.="<tr><td valign='top' style='background:#fff;color:#000;'>";
          $locationtype.=getTable($exon_rank,$frame,$snposincds,$snposinprotein,$blossumScore,$cai,
                                  $fb6_codon,$balck6_aa,$b6oneletter,$b6_codon,
                                  $b6_aa_polarity,$b6_hydropathy,$b6hydropathy,$balck6_aa,$SNP_HOME2,
                                  $fsnp_codon,$snp_aa,$oneletter,$snp_codon,$snp_aa_polarity,$snp_hydropathy,$snphydropathy,
                                  $snp_aa);
          $locationtype.="</td></tr>";
       }
       else{ 
           $locationtype.="<tr><td></td><td>$balck6_aa</td><td class='sectionOptions'>$snp_aa</td></tr>";}
    $locationtype.="</table>";
    return $locationtype ;         
}
/***************************************************************
   Given the user selection,
   This function generates the list of SNP 
   within the specified range
 displayRows($snplist,$rowcount,$page,$datatype,$selectedStrains,
                     $is_imputed,$selectedSources,$strainIndex,$snp_location,$fields){
  $snpcount=$snplist[0];$snpids=$snplist[1];
  global $MGI_GENEDETAIL;global $SNP_SEARCH_HOME;global $CGD_GBROWSE_37;global $SNPLOC; 
  $chr=$datatype{"chr"};$startpos=$datatype{"startpos"};$endpos=$datatype{"endpos"};
  $gene_id= $datatype{"geneid"};global $GENOTYPEMAP;
  $transcript_id= $datatype{"transcriptid"}; $item_type_id=$datatype{"qtype"};
  $i=1; $sourcelist=split(":",$selectedSources); $sourcecount=count($sourcelist);
****************************************************************/
function displayRowsCSV($snplist,$datatype,$selectedStrains,$is_imputed,$selectedSources,$strainIndex,$snp_location,$fields){
  $snpcount=$snplist[0];$snpids=$snplist[1];global $SNPLOC; global $MGI_GENEDETAIL;global $SNP_SEARCH_HOME;
  $chr=$datatype{"chr"};$startpos=$datatype{"startpos"}; global $CGD_GBROWSE_37;
  $endpos=$datatype{"endpos"};$gene_id= $datatype{"geneid"};
  $transcript_id= $datatype{"transcriptid"}; $item_type_id=$datatype{"qtype"};
  $i=0; $sourcelist=split(":",$selectedSources); $sourcecount=count($sourcelist);
  $selectedStrainlist=split(":",$selectedStrains); /**********/
  while(($i<$snpcount)){ // process every snp returned from user selection
         $snpid=trim($snpids[$i]);
         $all=0;$snpaccession_id=getSNPaccession($snpid,$selectedSources,$all);
         $snpdata=getsnpdata($con,$snpid);
         if(!empty($snpdata)){
            $chr=$snpdata{"chr"};$bp_position=$snpdata{"pos"};
            $black6_allele=$snpdata{"b6allele"};$snp_allele=$snpdata{"snpallele"};
            $is_CpG_site=$snpdata{"iscpg"};$is_intergenic=$snpdata{"istergenic"};
            $cpg="no"; if($is_CpG_site)$cpg="yes";
            $mutation="Transversion";$mutation_type=$snpdata{"mutation_type"};
            if($mutation_type==1)$mutation="Transition";$impgenotypedata=array();$genotypeData=array();
            $mgigene="N/A";$ensemt="N/A";$snplocation="Intergenic";$is_same=0;$thisgeno=""; global $GENOTYPEMAP;
            if($is_intergenic==0)$snplocation="";
            $strains=split(":",$selectedStrains);
           if(is_array($strains)&&(!empty($strains[0]))){
              if($is_imputed==1)$impgenotypedata=getGenotypedata($snpid,$selectedStrains,$selectedSources,$is_imputed);
              else{$genotypeData =getGenotypedata($snpid,$selectedStrains,$selectedSources,0);}
              $same_geno=1;$strain_count=count($selectedStrainlist);
              if(!empty($genotypeData))$same_geno=$genotypeData{"samegeno"};
              else $same_geno=$impgenotypedata{"samegeno"};
              $conf="H";$char="N";$html=0;$k=0; $l=3; // index of first strain
              while($k<count($selectedStrainlist)){
                   $strain_id=$strainIndex[$l]; ++$l;  ++$k;$conf="H";$char="N";//get this strain geno
                   if(!empty($impgenotypedata{"$strain_id"})){ //get the imputed genotype
                       $thisStraingeno=split(":",$impgenotypedata{"$strain_id"}); 
                       list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                       $thisgeno.=",$geno";
                   }
                   else if(!empty($genotypeData{"$strain_id"})){ // for other sources
                       $thisStraingeno=split(":",$genotypeData{"$strain_id"});
                       list($has_conflict,$geno)=split(":",getThisGenotype($thisStraingeno,$html));
                       $thisgeno.=",$geno";
                  }
                  else{// this strain does not have a genotype data /strain not provided by this source
                      $conf="H";$char="N";
                      if($strain_id==7){$thisgeno .=",$black6_allele";}
                      else{$thisgeno .=",N";} 
                   }
               }
            } // end of if(is_array($strains)&&(!empty($strains[0]))){
            // get minor allele frequency :
           
            $impgenoSum="";$genoSumarry=""; $minorallele=""; $genomap="";
            if($sourcecount==1)
               $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);
             else{
               $genoSumarry= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,0);
               if($is_imputed==1)$impgenoSum= getMinorAlleleFreq($snpid,$black6_allele,$selectedSources,$is_imputed);
               
             }
            $genomap= parseMinorAllele($genoSumarry,$impgenoSum,$snp_allele,$black6_allele);
           /*
              now process the genomap to compute the minor allele frequency range.
              MAF= (minAllele_count/(totalStrain_count-missingGenotype_count))*100
              $minor_allele_strain_freq =ceil(($minor_allele_strain_freq/($NoNs_count-$Ns))*100);
           */
            $minallelef=computeMAF($genomap,$black6_allele,$snp_allele,$sourcelist);
           $minal  =  $minallelef{"minal"};$minmaf =  $minallelef{"minmaf"};
           $nmap   =  $minallelef{"nmap"};$totalStrain = $minallelef{"totalStrain"};$minorallele="$minal - $minmaf "; 
           if($is_intergenic==0){$snplocation="";
               $snplocdata=getSNPlocation($snpid,$gene_id,$transcript_id,$snp_location);
               while(list($gene_id,$transcriptdata)=each($snplocdata)){
                     $transcriptdatal=split(":",$transcriptdata);
                     $gbrowse_start=$bp_position -200; $gbrowse_end=$bp_position+200;                               
                     $row="Chr$chr,$bp_position";$row.=$thisgeno; $loccount=count($transcriptdatal);
                     $row.=",$snpaccession_id";
                     if($fields{"v_type"}==1){$row.=",$mutation";}
                     if($fields{"g_name"}==1){
                        $mgigene=getMgiGene($gene_id,0);
                        $mgigene_symbol=$mgigene[1];$mgigene_id=$mgigene[0];
                        if(empty($mgigene))$mgigene_symbol=getEnsGene($gene_id);
                        $row.=",$mgigene_symbol";
                     }$o=0;
                     if($fields{"f_class"}==1){
                        while($o<($loccount-1)){
                                $trans_location=$SNPLOC{$transcriptdatal[$o]};
                                if($transcriptdatal[$o]==1)$trans_location="CSyn";
                                else if($transcriptdatal[$o]==2)$trans_location="CNSyn";  
                                else if($transcriptdatal[$o]==5)$trans_location="StopGained";   
                                else if($transcriptdatal[$o]==6)$trans_location="StopLost";
                                else if($transcriptdatal[$o]==7)$trans_location="StartGained"; 
                                else if($transcriptdatal[$o]==0)$trans_location="Intronic";
                                else if($transcriptdatal[$o]==3)$trans_location="3'UTR";
                                else if($transcriptdatal[$o]==4)$trans_location="5'UTR";
                                else if($transcriptdatal[$o]==-3)$trans_location="Exon(non protein-coding)";
                                else if($transcriptdatal[$o]==-2)$trans_location="Intron(non protein-coding)";
                                else if($transcriptdatal[$o]==8)
                                      $trans_location="StartLost";                                        
                                if($o==0)$snplocation ="$trans_location";
                                else{$snplocation .="-$trans_location";} 
                                ++$o;                          
                           
                          }$row.=",$snplocation";
                      }
                     if($fields{"cpg"}==1){$row.=",$cpg";}
                      $row.=",$black6_allele/$snp_allele";
                     if($fields{"m_allele"}==1){$row.=",$minorallele";}
                     if($fields{"m_genotype"}==1){$row .=",$nmap";}
                      $row.="\015\012";
                      print"$row";
               } //end of  while(list($gene_id,$transcriptdata)=each($snplocdata)){
          } // end of  if($is_intergenic==0){
          else{
               $row="Chr$chr,$bp_position";
               $row.=$thisgeno;$row.=",$snpaccession_id";
               if($fields{"v_type"}==1){$row.=",$mutation";}
               if($fields{"g_name"}==1){$row.=",$mgigene";}
               if($fields{"f_class"}==1){$row.=",$snplocation";}
               if($fields{"cpg"}==1){$row.=",$cpg";}
               $row.=",$black6_allele/$snp_allele";
               if($fields{"m_allele"}==1){$row.=",$minorallele";}
               if($fields{"m_genotype"}==1){$row .=",$nmap";}
              print"$row\015\012";
          }
      } //end of if snpdata
      ++$i;    
   }
}

?>
