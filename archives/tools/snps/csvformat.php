<?php
 
  $filename="query_result.csv";
  // header("Content-Disposition: attachment; filename=\"$filename\"");
  header("Content-type: application/octet-stream");
  header("Content-Disposition: attachment; filename=\"$filename\"");
  header("Cache-Control: cache, must-revalidate");
  header("Pragma: public");
  header("Content-type: text/plain; charset=UTF-8");

  $newpath= ini_get('include_path');
  $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils";
  $newpath .=":/raid/Users/LucieH/utils";
  
  ini_set('include_path', $newpath);
  require_once("build37/global_utils.php");
  require_once("build37/cgd_utils.php");
  global $SNPLOC; global $CGDTYPE;
  if(empty($SNPLOC))loadSNPloc();
 
  global $SNPLOC; global $CGDTYPE;
  if(empty($SNPLOC))loadSNPloc();
 
        $mutation_type=0;$cpgsite_check="";$confScore="";global $outype;$rowcount=0;$page=0;
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
        $page=filterInput(getSelectedItemList($myquerystring,"page"));$querystring="";
        $snp_loc=filterInput(getSelectedItemList($myquerystring,'snp_loc'));
        $genotypeChek=filterInput(getSelectedItemList($myquerystring,"genotype"));
        $qType=filterInput(getSelectedItemList($myquerystring,'qType'));
        $transcript_local_id=0;$transcript_local_id=filterInput(getSelectedItemList($myquerystring,"tid"));
        $datatype=ValidateData($singlequery);
  $error=$datatype{"error"};
 if(!empty($error)){
    print" $error \015\012 You Can query using the following terms:\015\012";
    print"Single chromosome:bpPosition (bp or Mbp) :4:20266796 OR 4:20.266796\015\012
          Genomic range, example : chromosome:Position1-Position2 (bp or Mbp) : 4:20231103-20575856 OR 4:20.2-20.5\015\012
          SNP accession id :</b>NES08626062,GNF3-103994 ,rs6356141\015\012
          Ensembl Ids :  gene_id, transcript, protein_id\015\012
          MGI genes  : gene_id, gene_name (Zfy1,..) \015\012
          Entrez Gene : gene_id, gene_name (LOC100034724, ...)\015\012";
 }
  else{
      $querystring=$_SERVER['QUERY_STRING'];
      $myquerystring= split("&",$_SERVER['QUERY_STRING']);
      $selectedStrains=getSelectedItemList($myquerystring,"st");
      $selectedSources= getSelectedItemList($myquerystring,"source");
      if(empty($selectedSources))$selectedSources="21"; //dfault sources (Perlegen and Imputed)
      if($datatype{"qtype"}==11){ #query by snp accession id, add snp source to the list
         $sc_id=$datatype{"qsource"};
         $selectedSources.=":$sc_id";
      }
      $sourcedata=getSNPsource($selectedSources);
      $sourcelist=$sourcedata{"sourcelist"};
      $is_imputed= $sourcedata{"imputedflag"};
      processQueryCSV($datatype,$snp_loc,$genotypeChek,$selectedSources,$confScore,$selectedStrains,$rowcount,$page,$fields,$querystring);      
   }
  
?>
  
