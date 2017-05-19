<?php

   $newpath= ini_get('include_path');
   $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils";
   $newpath .=":/raid/Users/LucieH/utils";
   
   ini_set('include_path', $newpath);
   require_once("build37/global_utils.php");
   $source_id=21; //default imputed source
   $querystring=$_SERVER['QUERY_STRING'];   
  
   ///////////////////////
   $mutation_type=0;$cpgsite_check="";$confScore="";global $outype;$rowcount=0;$page=0;$filterlink="";
   $fields{"v_type"}=0;$fields{"g_name"}=0;$fields{"f_class"}=0;$fields{"cpg"}=0;$fields{"sg"}=0;
   $myquerystring= split("&",$_SERVER['QUERY_STRING']);
   if(empty($singlequery))$singlequery=filterInput(getSelectedItemList($myquerystring,"query_single"));
    $fields{"sg"}=filterInput(getSelectedItemList($myquerystring,"sg"));
    $fields{"m_allele"}=0;$fields{"m_genotype"}=0;
    $fields{"v_type"}=filterInput(getSelectedItemList($myquerystring,"vt"));
    $fields{"f_class"}=filterInput(getSelectedItemList($myquerystring,"fc"));
    $fields{"cpg"}=filterInput(getSelectedItemList($myquerystring,"cp"));
    $fields{"g_name"}=filterInput(getSelectedItemList($myquerystring,"gn"));
    $fields{"m_allele"}=filterInput(getSelectedItemList($myquerystring,"ma"));
    $fields{"m_genotype"}=filterInput(getSelectedItemList($myquerystring,"mg"));
    $confScore =filterInput(getSelectedItemList($myquerystring,"confScore"));
    $geneid=0; $geneid=filterInput(getSelectedItemList($myquerystring,"gid"));
    $rowcount=filterInput(getSelectedItemList($myquerystring,"rcount"));
    $page=filterInput(getSelectedItemList($myquerystring,"page"));
    $snp_loc=filterInput(getSelectedItemList($myquerystring,'snp_loc'));
    $genotypeChek=filterInput(getSelectedItemList($myquerystring,"genotype"));
    $qType=filterInput(getSelectedItemList($myquerystring,'qType'));
    $transcript_local_id=0;$transcript_local_id=filterInput(getSelectedItemList($myquerystring,"tid"));
    $datatype=ValidateData($singlequery);
    $error=$datatype{"error"};
//////////////////////////////////////////////////////////////
   $result="<a href='$SNP_SEARCH_HOME?$querystring' id='hq'>hd</a>";
   $test="<a href='$SNP_SEARCH_HOME?$querystring' id='ta'>hd:$SNP_SEARCH_HOME - $querystring</a>";
   if(!empty($error)){
       $result.='<table id="maintable"><tr>
                       <td id="leftp" valign="top" width="10"></td><td id="middlep" valign="top">';
       $result.=" <div id='qhint' ><span>$error </span><br/><h2>You Can query using the following terms:</h2>";
       $result.='      <ul> <li> <b>Single chromosome:bpPosition (bp or Mbp) :</b>4:20266796 OR 4:20.266796</li>
                                   <li> <b>Genomic range, example : chromosome:Position1-Position2 (bp or Mbp) : </b>
                                       4:20231103-20575856 OR 4:20.2-20.5</li>
                                   <li><b>SNP accession id :</b>NES08626062,GNF3-103994 ,rs6356141</li>
                                   <li><b>Ensembl Ids : </b> gene_id, transcript, protein_id</li>
                                   <li><b>MGI genes  : </b> gene_id, gene_name (Zfy1,..) </li>
                                  <li><b>Entrez Gene  :</b> gene_id, gene_name (LOC100034724, ...)</li>
                           </ul></div></td></tr></table>';
   }
   else{
       ///////////////////////////////////
       $selectedStrains=getSelectedItemList($myquerystring,"st");
       if(empty($selectedStrains))$genotypeChek=0;
       $strainsList=getSelectedStrainList($selectedStrains);
       $selectedSources= getSelectedItemList($myquerystring,"source");
       if(empty($selectedSources))$selectedSources="1:21"; //dfault sources (Perlegen and Imputed)
       if($datatype{"qtype"}==11){                         #query by snp accession id, add snp source to the list
          $sc_id=$datatype{"qsource"};
          $selectedSources.=":$sc_id";
        }
       $sourcedata=getSNPsource($selectedSources);
       $sourcelist=$sourcedata{"sourcelist"};
       $sg=$fields{"sg"};
       $filterlink="$SNP_SEARCH_HOME?query_single=$singlequery&amp;sg=$sg&amp;genotype=$genotypeChek&amp;confScore=$confScore";
       $filterlink .="&amp;rcount=$rowcount";
       $sources=split(":",$selectedSources);
       foreach($sources as $source_id){$source_id=trim($source_id);$filterlink.="&amp;source=$source_id";}
       $strains=split(":",$selectedStrains);
       if(is_array($strains)&&!empty($strains[0]))
          foreach($strains as $strain_id)$filterlink .="&amp;st=$strain_id";
      
       $is_imputed= $sourcedata{"imputedflag"};$chr=$datatype{"chr"};$startpos=$datatype{"startpos"};
       $endpos=$datatype{"endpos"};
       $gene_id= $datatype{"geneid"};
       $query_single=$datatype{"query_single"};
       $transcript_id= $datatype{"transcriptid"}; 
       $item_type_id=$datatype{"qtype"};
       $qtype_name=$CGDTYPE{$item_type_id}{"type_name"};
       $chrStr=$chr;if($chr==20)$chrStr="X"; if($chr==21)$chrStr="Y"; if($chr==22)$chrStr="Mt";
       
       $snplocationdata=setSNPloc($snp_loc); // get the user specified genomic region and annotations
       $locCount=$snplocationdata[0];$snp_location=$snplocationdata[1];$location_name=$snplocationdata[2];
       $is_intergenic=$snplocationdata[3];
       
       $SNPlist= getSNPsIDs($con,$chr,$startpos,$endpos,$selectedSources,$snp_location,
                       $is_intergenic,$gene_id,$transcript_local_id,$locCount);
       $currenttotalSNP=$SNPlist[0]; $type=0;
       $result.=displayResulttable($currenttotalSNP,$strainsList,$SNPlist,$rowcount,$page,$datatype,$selectedStrains,
                        $is_imputed,$selectedSources,$strainIndex,$snp_loc,$fields,$querystring,$filterlink,$type);
   }
   echo $result;
   //echo "<div>--->$filterlink --></div>";
  
?>
