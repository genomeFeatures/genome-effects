<?php  //redirect to csvformat if the output type is not htm  
   session_start(); 
   error_reporting(0); //disable error messages 
    if(!empty($_GET["outype"])) $outype =$_GET["outype"];
    $querystring=$_SERVER['QUERY_STRING'];
    if($outype>1)
       header("Location:http://cgd-dev.jax.org/cgdsnpdb/csvformat.php?$querystring");
 ?>
<?php 
  $newpath= ini_get('include_path');
  $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils";
  $newpath .=":/raid/Users/LucieH/utils";
 
  ini_set('include_path', $newpath);
  require_once("build37/global_utils.php");
  require_once("build37/cgd_utils.php");

  global $SNPLOC; global $CGDTYPE;global $ASSEMBLY_BUILD;global $DBVERSION;
  if(empty($SNPLOC))loadSNPloc();global $STRAINGROUPS;
  if(empty($STRAINGROUPS))loadStrainGroups(); 
  $title="CGD SNP Database $DBVERSION ($ASSEMBLY_BUILD) -- The Jackson Lab";
 ?> 
 <!DOCTYPE html
 PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
     "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title><?php  print"$title";?></title>
<!-- jax header -->
 <meta http-equiv="Content-Language" content="en-us"/>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>
  <meta name="Classification" content="research"/>
  <meta name="verify-v1" content="Asn8D9XHy175bBjQrJtq2A5Yx34ANQJaiVUpdHkxmEI=" />
  <meta name="Description" content="CGDSNP is a high quality Mouse SNP database with more 
           than 8 Millions SNPs from  74 strains of laboratory mice."/>
  <meta name="Keywords" content="SNPs,mouse SNPs, SNP database,SNP,mouse, SNP dataset,SNP data"/>
<?php global $SNP_HOME; global $CGD_HOME;
  $links="<link rel=\"stylesheet\" type=\"text/css\" href=\"$SNP_HOME/css/cgd.css\" />";
  $links.="<link rel=\"shortcut icon\" href=\"$CGD_HOME/favicon.ico\" />";
  $links.="<script type=\"text/javascript\" src=\"$SNP_HOME/js/cgdjavascript.js\"></script>";
  $links.="<link rel=\"stylesheet\" type=\"text/css\" href=\"$SNP_HOME/css/movecolumns.css\" />";
  $links.="<script type=\"text/javascript\" src=\"$SNP_HOME/js/cgdjavascript.js\"></script>";
  $links.="<script type=\"text/javascript\" src=\"$SNP_HOME/js/movecolumns.js\"></script>";
  echo $links;
 ?>
<script type="text/javascript">
 <!--
   function showThisItem(itemid,displayType){
       document.getElementById(itemid).style.display=displayType;
   }
   function hideThisItem(itemid){
      document.getElementById(itemid).style.display="none";
   }
   //--> 
 </script>
<style type="text/css">
  #qhint{width:90%; color:#33A2A2; font-size:1.1em;}
  #qhint ul{list-style: circle; margin-left:2.5em;}
  #qhint h2{color:#3F2C6D;}
  #qhint b{color:#333;}
  #qhint span{padding:2em; font-size:0.9em;color:red;}
  fieldset{margin:0;padding:0; 
    -moz-border-radius: 5px;    /* For Mozilla Firefox */
    -khtml-border-radius: 5px;  /* For Konqueror */
    -webkit-border-radius: 5px; /* For Safari */
    border-radius: 5px;         /* For future native implementations */
    background:#eee;
    border:solid 1px #eee;
  }
  fieldset .fh{width:100%;margin:0;font-weight:bold;taxt-align:center;background:#3F2C6D;color:#eee; font-size:1.3em;}
  fieldset p,.tdata,.sectionName,.strainlabel{font-size:0.85em;}
  .sr{padding:0.0em; margin:0em; height:auto; overflow:auto; width:auto;}
  .sr h2{color:#3F2C6D; margin:0; padding:0;font-size:1.1em;white-space: nowrap;}
  .sr span{font-size:0.95em; margin:0;}
  .ah,.am,.al,.th,.tm,.tl,.ch,.cm,.cl,
  .gh,.gm,.gl,.h,.n{margin:0;padding:0.1e; font-size:0.85em;}
  .strainlabel{background:rgb(218,233,218);}
  .ah,.am,.al,.th,.tm,.tl,.ch,.cm,.cl,.gh,.gm,.gl{color:#fff;}
  /***** Yellows *****************/
  .gh{background:rgb(215,215,0); } 
  .gm{background:rgb(215,215,150);} 
  .gl{background:rgb(215,215,200);}
  
  /***** Blues ********************/
  .ch{background:rgb(0,51,204);} 
  .cm{background:rgb(153,173,235);}
  .cl{background:rgb(204,214,245);}
  
  /***** Black ********************/
  .ah{background:rgb(51,51,51);} 
  .am{background:rgb(102,102,102);} 
  .al,.sectionName{background:rgb(204, 204, 204);}
  
  /***** Red ********************/
  .th{background:rgb(255,0,0);} 
  .tm{background:rgb(255,128,128);} 
  .tl{background:rgb(255,204,204);}
  
   .h,.n,{background:#fff;}
   .h{rgb(0,51,204);}
   .n,.strainlabel{color:#333;}
   .strainlabel{padding-left:0.5em;padding-right:0.5em;padding-bottom:0.6em;}
   .strainlabel ul{margin:0;padding:0;list-style:none;}
   .strainlabel li{margin:0;padding:0;float:left;clear:both;height:0.9em;}
   #snptable {margin-left:1em;}
   #snptable td{padding-left:0.6em;}
  .sr ul{margin:0; padding:0;margin-top:0.4em;list-style:none; font-size:0.9em;}
  .sr li{white-space: nowrap ;  padding-left:0.3em;}
  .sectionOptions{padding:0.4em;font-size:0.9em;background:rgb(218,233,218);}
  .dtable{width:90%; margin:auto;}
  .dtable .sectionName{background:rgb(212,212,212);}
  .abutton{font-size:0.7em;float:right;}
  #hfield td{padding-left:1.2em;margin:0;font-size:0.75em;}
  .check{margin:0;padding:0; font-size:0.4em;}
  #tresult table{background:#fff;}
  #hq{display:none;}
 </style>
  <script type="text/javascript">
    <!--
   function show_hide_column(col_no, do_show) {
    var stl;
    if (do_show) stl = 'block'
    else {
        stl = 'none';
        // add this column to the hidden fields list
       // document.getElementById("hfield")
     }

    var tbl  = document.getElementById('snptable');
    var rows = tbl.getElementsByTagName('tr');
  
    for (var row=0; row<rows.length;row++) {
      var cels = rows[row].getElementsByTagName('td')
      cels[col_no].style.display=stl;
    }
    
  }

    //-->
    </script>
</head>
<body>
  <div class="mainWindow2">
    <form name='mainSearch' action='parseQuery.php?rtn=' method='get'>
     <?php
        $title="Center for Genome Dynamics Mouse SNP Database";
        displayPageHeader($title);
        displayTopNav();
        $filename="";$geneloChr="";$singlequery="";$singlequeryOption="";$snp_loc=12;$SNPsource="";                          
        $snpid=0;
        if(!empty($_GET["snpid"])){
            $snpid=filterInput($_GET["snpid"]);
            $selectedStrains=filterInput($_GET["ss"]);
            displaySnpDetailView($snpid,$selectedStrains);
        } 
       else{
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
        //if(isset($_SESSION['genome']))echo("Session set");
        if(!empty($error)){
            print'<table id="maintable"><tr>
                       <td id="leftp" valign="top" width="10"></td><td id="middlep" valign="top">';
            print" <div id='qhint' ><span>$error </span><br/><h2>You Can query using the following terms:</h2>";
            print'      <ul> <li> <b>Single chromosome:bpPosition (bp or Mbp) :</b>4:20266796 OR 4:20.266796</li>
                             <li> <b>Genomic range, example : chromosome:Position1-Position2 (bp or Mbp) : </b>
                                       4:20231103-20575856 OR 4:20.2-20.5</li>
                             <li><b>SNP accession id :</b>NES08626062,GNF3-103994 ,rs6356141</li>
                             <li><b>Ensembl Ids : </b> gene_id, transcript, protein_id</li>
                             <li><b>MGI genes  : </b> gene_id, gene_name (Zfy1,..) </li>
                             <li><b>Entrez Gene  :</b> gene_id, gene_name (LOC100034724, ...)</li>
                     </ul></div></td></tr></table>';
          }
         
          else{
            $querystring=$_SERVER['QUERY_STRING'];
            $myquerystring= split("&",$_SERVER['QUERY_STRING']);
            $selectedStrains=getSelectedItemList($myquerystring,"st");
            if(empty($selectedStrains))$genotypeChek=0;
            $selectedSources= getSelectedItemList($myquerystring,"source");
            if(empty($selectedSources))$selectedSources="21"; 
                                             //default sources (Perlegen and Imputed)
            if($datatype{"qtype"}==11){     #query by snp accession id, add snp source to the list
               $sc_id=$datatype{"qsource"};
               $selectedSources.=":$sc_id";
             }
            $sourcedata=getSNPsource($selectedSources);
            $sourcelist=$sourcedata{"sourcelist"};$is_imputed= $sourcedata{"imputedflag"};
            if($snpid>0){
               $selectedStrains=$_GET["ss"];
               displaySnpDetailView($snpid);
             }
            else{
               //$chr=$datatype{"chr"};$startpos=$datatype{"startpos"};$endpos=$datatype{"endpos"}; 
               processQuery($datatype,$snp_loc,$genotypeChek,$selectedSources,$confScore,
                        $selectedStrains,$rowcount,$page,$fields,$querystring); 
            }       
         }
       }
       ?>          
     <?php displayFooter(); ?>
    </form>
   </div> </body></html>
