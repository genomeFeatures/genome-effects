<?php
/*********************************************************** 
  This file stores 
  all global variables and functions to be used on the entire site
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : December 2009

***************************************************************/ 
 $newpath= ini_get('include_path');
  $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils"; 
  $newpath .=":/xraid/Users/LucieH/utils";
 
  ini_set('include_path', $newpath);
  require_once("build37/user.php");
  $Mbp =1000000;
  $RANGEMAX=10*$Mbp;
  $RANGEMAXs="10 Mbp";
  $SNP_SEARCH_HOME="http://cgd.jax.org/cgdsnpdb/search";
  $SNP_HOME="http://cgd.jax.org/cgdsnpdb";
  $SNP_HOME2="http://cgd.jax.org/cgdsnpdb";
  $CGD_HOME="http://cgd.jax.org";
  /*****************************
    $SNP_SEARCH_HOME="http://demon.jax.org/lnh/prototype/lnh_cgd/cgdsnpdb_wi/search";
    $SNP_HOME="http://demon.jax.org/lnh/prototype/lnh_cgd/cgdsnpdb_wi";
  ***************************************/
  $ENSEMBL_46_36_PAGE="http://aug2007.archive.ensembl.org/Mus_musculus/exonview?db=core;transcript=";
  $ENSEMBL_47_37_PAGE="http://www.ensembl.org/Mus_musculus/Transcript/Summary?db=core;t=";
  $ENSEMBL_47_37_PAGE="http://www.ensembl.org/Mus_musculus/Gene/Variation_Gene/Table?db=core;g=";
  //$ENSEMBL_47_37_PAGE="http://oct2007.archive.ensembl.org/Mus_musculus/geneview?gene=";
  $ENSEMBL_47_37_TRANSCRIPT="http://www.ensembl.org/Mus_musculus/Transcript/Summary?db=core;t=";
  $ENSEMBL_47_37_PEPTIDE="http://www.ensembl.org/Mus_musculus/Transcript/ProteinSummary?db=core;t=";
  $ENSEMBL_CURRENT_PAGE="http://www.ensembl.org/Mus_musculus/Gene/Variation_Gene/Table?db=core;g=";
  $ENSEMBL_CURRENT_TRANSCRIPT="http://www.ensembl.org/Mus_musculus/Transcript/Summary?db=core;t=";
  $ENSEMBL_CURRENT_PEPTIDE="http://www.ensembl.org/Mus_musculus/Transcript/ProteinSummary?db=core;t=";
  $CGD_GBROWSE ="http://cgd.jax.org/cgi-bin/gbrowse/mouse";
  $CGD_GBROWSE_37="http://cgd.jax.org/cgi-bin/gbrowse/mouse37";
  $MGI_GENEDETAIL="http://www.informatics.jax.org/searches";
  $ENTREZ_GENE="http://www.ncbi.nlm.nih.gov/sites/entrez?db=gene&term=";
  $BUILD36_SNP_HOME="http://cgd.jax.org/cgdsnpdb/build_36";
  $BUILD36_IMPUTED_DOWNLOAD="http://cgd.jax.org/datasets/popgen/imputed36.shtml";
  $BUILD37_IMPUTED_DOWNLOAD="http://cgd.jax.org/datasets/popgen/imputed.shtml";
  $NCBI_SNP_PAGE="http://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=";
  $ASSEMBLY_BUILD="NCBI Mouse Build 37";
  $DBVERSION ="v1.5";
  $SNPCODON[0]=0;
  $UNIT="bp";
  $SNPLIST; //array of SNPs in a given genomic region

$SNP_SEARCH_HOME ="http://cgd.jax.org/cgdsnpdb/parseQuery.php";
#$SNP_SEARCH_HOME ="http://cgd-dev.jax.org/cgdsnpdb/parseQuery.php";

$SNP_SEARCH_CSV="http://cgd.jax.org/cgdsnpdb/csvformat.php";
#$SNP_SEARCH_CSV="http://cgd-dev.jax.org/cgdsnpdb/csvformat.php";
$STRAINSBYSOURCE;  //stores the mapping between a
                  // given source and the corresponding strains count
              
/***********************************************************  
   stores the mapping between a
   given source and the corresponding SNP count 
************************************************************/
$SNPSBYSOURCE{1} ="8,224,916";     //perlegen SNP COUNT 
$SNPSBYSOURCE{20} ="549,645";     //MusDiv Diversity Array - Yang et al.2011
$SNPSBYSOURCE{21} ="65,243,635";  //Imputed - Jeremy R. Wang et al. 2012
$SNPSBYSOURCE{16} ="7,867,995";    //Imputed SNP COUNT
$STRAINGROUPS;
/*$SNPSBYSOURCE{2} ="138,594";       //BROAD snp count
//$SNPSBYSOURCE{6} ="155,749";       //GNF
$SNPSBYSOURCE{13} ="2,122,059";    //Celera
$SNPSBYSOURCE{16} ="7,867,856";    //Imputed SNP COUNT 
$SNPSBYSOURCE{17} ="548,363";      //Mus Diversity Array snp count
$SNPSBYSOURCE{18} ="24,608";       //Paigen
$SNPSBYSOURCE{19}="667";           //Wild derived

/***********************************************************  
   stores the mapping between a
   given source and the corresponding source name 
************************************************************/
$SNPSOURCE{21} ="Imputed - Jeremy R. Wang et al. 2012";             //Imputed 
$SNPSOURCE{20} ="MusDiv Diversity Array - Yang et al.2011";       //Mus Diversity Array snp count
$SNPSOURCE{1} = "NIEHS";             //perlegen SNP COUNT 
$SNPSOURCE{16} ="Imputed - Szatkiewicz et al.2008";             //Imputed

/*$SNPSOURCE{6} = "GNF";              //GNF
$SNPSOURCE{18} ="Paigen";          //Paigen
$SNPSOURCE{13} ="Celera";         //Celera
$SNPSOURCE{19} ="Wild Derived";   //Wild derived
$SNPSOURCE{2} ="Broad";           //BROAD snp count
*/
/********************* CC founders strain map ******************/
$CCFOUNDER{1}="129S1/SvlmJ";
$CCFOUNDER{2}="A/J";
$CCFOUNDER{7}="C57BL/6J"; 
$CCFOUNDER{8}="CAST/EiJ";
//$CCFOUNDER{13}="NOD/LtJ";
$CCFOUNDER{129}="NOD/ShiLtJ"; 
$CCFOUNDER{16}="WSB/EiJ";                                                                                                   $CCFOUNDER{213}="NZO/HlLtJ";
$CCFOUNDER{214}="PWK/PhJ";

$CGDTYPE;
$ENSEMBLGENES;
$ENSEMBLTRANSCRIPTS;
$ENSEMBLPROTEINS;
$SNPLOC;
$SNPLIST;
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
/*************************************************************
   stores the cgd type mapping
**************************************************************/
function loadCgdType(){
  global $CGDTYPE;
  global $con;
  if(!$con)$con=getConnection();
  $query="select items_type_id,type_name,item_table_name,item_key_name from cgd_items_type";
  $query .=" where items_type_id>6";
  $result = mysql_query($query,$con);
  $numRows = mysql_numrows($result);
  while($row=mysql_fetch_assoc($result)){
        $items_type_id= $row["items_type_id"];
        $type_name = $row["type_name"];
        $item_table_name= $row["item_table_name"];
        $item_key_name = $row["item_key_name"];
        $CGDTYPE{$items_type_id}{"type_name"}=$type_name;
        $CGDTYPE{$items_type_id}{"item_table_name"}=$item_table_name;
        $CGDTYPE{$items_type_id}{"item_key_name"}=$item_key_name;
  }
  $CGDTYPE{0}{"type_name"}="Genomic location";
  $CGDTYPE{0}{"item_table_name"}="";
  $CGDTYPE{0}{"item_key_name"}="";      
}
/*************************************************************
   stores SNP location mapping
**************************************************************/
function loadSNPloc(){
  global $SNPLOC;
  global $con;
  if(!$con)$con=getConnection();
  $query="select _loc_func_key,description  from snp_loc_func";
  $result = mysql_query($query,$con);
  $numRows = mysql_numrows($result);
  while($row=mysql_fetch_assoc($result)){
        $_loc_func_key= $row["_loc_func_key"];
        $description = $row["description"];
        $SNPLOC{$_loc_func_key}= $description;
  }  
  $SNPLOC{"10"}= "SNPs on Exons"; 
  $SNPLOC{"12"}= "All SNPs";        
}


/**************************************************************
  displays the footer menue 
***************************************************************/
function displayFooter(){
 print'<div id="footer">';
   displayTopNav();
 print'</div>';
} 
/**************************************************************
  displays contact navigation menue 
***************************************************************/
function displayTopNav(){
  global $CGD_HOME; global $SNP_HOME; global $BUILD36_SNP_HOME;
 //$BUILD36_SNP_HOME="http://cgd.jax.org/cgdsnpdb/build_36";
  print"<div class='topnav'><center><ul><li><a href='$SNP_HOME/'>";
  print'CGDSNP Home </a></li><li class="lborder_w">';
  print"<a href='$CGD_HOME/' target='_blank'>Center HOME</a> </li>";
  print"<li class='lborder_w'><a href=\"$SNP_HOME/utils/contact.php\"> Contact Us</a> </li>";
  print "<li class='lborder_w'><a href=\"$SNP_HOME/utils/collaborators.php\">Collaborators</a></li>";
  print" <li class='lborder_w' style='width:auto;' onmouseover=\"showThisItem('arc','block')\" onmouseout=\"hideThisItem('arc');\">
      <a href=\"$BUILD36_SNP_HOME\" >Archives <span style='font-size:0.75em;'>(NCBI Build 36)</span></a>";
  print"</li></ul></center></div>";
}

/****************************************************************
  displays the page header section of each generated html file
*****************************************************************/
function displayPageHeader($title)
{
 global $SNP_HOME;
 print'<div id="pheader">
        <table><tr><td id="logo"><a href="http://www.jax.org" target="_blank">';
 print"         <img src='$SNP_HOME/images/CGD-web-logo.gif' class='cgdlogo' alt='Jax logo'/> </a></td>
        <td><div id='pagetitle'><h1>$title</h1></div></td>
             </tr></table>
     </div>";
 
}

/***************************************************************
  displays the browse navigation menue 
 ***************************************************************/
function displayLeftNav(){
  global $SNP_HOME;
  print'<fieldset id="news"><div class="fheader"> Browse</div>';
  print'   <ul style="padding-left:0.5em;">';
  print"     <li><a href=\"$SNP_HOME\" >CGDSNP Home</a></li>";
  print"     <li><a href=\"$SNP_HOME/build37/data/imputedsnps\" target='_blank'>Downloads</a></li>";
  print"     <li><a href=\"$SNP_HOME/utils/snp_data_report.php\">About the Database</a></li>";
  print"     <li><a href=\"$SNP_HOME/utils/public_mysql.php\">public MySQL Server</a></li> ";
  print"     <li><a href=\"$SNP_HOME/utils/snp_strain_by_source.php\">Strains By Source</a></li>";
  print"     <li><a href=\"$SNP_HOME/docs/user_doc.php\">How to use this tool</a></li></ul>";
  print"  </fieldset>";
}
function displaySummary(){
  print' <fieldset id="itable"><div  class="fheader">Build 37 SNPs Summary</div>
                  <table style="padding-left:0.5em;">
                     <tr><td><b>Total : &nbsp;&nbsp;</b>
                     <span class="sectionOptions">66,028,809</span></td></tr>
                     <tr><td><b>Transition:&nbsp;&nbsp;</b>
                             <span  class="sectionOptions">42,744,395(65%)</span>
                     </td></tr>
                    <tr><td><b>Transversion:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">23,284,436(35%)</span>
                    </td></tr>
                     <tr><td><b>Intergenic:&nbsp;&nbsp;</b>
                           <span  class="sectionOptions">38,919,578(59%) </span>
                     </td></tr>
                     <tr><td><b>On Transcript:&nbsp;&nbsp;</b>
                          <span class="sectionOptions">27,109,231(41%)</span>
                     </td></tr>
                     <tr><td><b>Intronic:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">26,102,091</span>
                     ; <b>Exonic:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">1,102,004</span>
                     </td></tr>
                     <tr><td><b>UTR:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">699,677</span>
                     ;<b>Synonymous:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">290,014</span></td></tr>
                     <tr><td><b>Non-Synonymous:&nbsp;&nbsp;</b>
                          <span  class="sectionOptions">161,260</span></td></tr>
          </table> </fieldset>';
}
/******************************************************************
  displays the archive section 
*******************************************************************/
function displayArchive(){
 global $BUILD36_SNP_HOME;
 print'<fieldset id="arc"><div class="fheader">Archives</div>';
 print"<ul><li><a href=\"$BUILD36_SNP_HOME\">NCBI Build 36 SNPs</a></li></ul>
       </fieldset>";
}
/******************************************************************
 displays the announcement section
 ******************************************************************/
function displayAnnouncement(){
 print' <fieldset id="announce">
      <div class="fheader">Announcements</div>
      <ul id="mlist">
          <li><h3>August 2012</h3>
              <ul><li> SNP database v1.5 maintenance release</li>
              </ul>
          </li>
       </ul>
   </fieldset>';
}
/*************************************************************** 
   This function generates the strains list mapping
****************************************************************/   
function getStrainMainList($con)
{
        global $con;
        if(!$con)$con=getConnection();
        $query="select distinct s.strain_id,strain_name from snp_strain s,";
        $query .=" snp_strain_by_source sc where s.strain_id =sc.strain_id";
        $query.="  order by strain_name";
        $result = mysql_query($query,$con);
        $numRows = mysql_numrows($result);
        $data="";
        while($row=mysql_fetch_assoc($result)){
              $strainID= $row["strain_id"];
              $strainName = $row["strain_name"];
              $data{$strainID}=$strainName;
          }        
       return $data;
}
/*************************************************************** 
   This function generates the sources list 
    for a given strain or returns all the sources if strainid<=0
****************************************************************/
function getTheSourceList($strainid){
   global $con;$sources=""; global $SNPSOURCE;
   //if($strainid<=0) return $SNPSOURCE;
   if(!$con)$con=getConnection();
   $query="select distinct s.source_id,s.source_name from snp_source s,snp_strain_by_source sc  ";
   $query.=" where sc.source_id=s.source_id";
   if($strainid>0)$query .=" and sc.strain_id =$strainid ";
   $query .=" order by source_name";
   $result = mysql_query($query,$con) or die("bad query :$query");
   $numRows = mysql_numrows($result);
   while($row=mysql_fetch_assoc($result)){
      $source_name= $row["source_name"];
      $source_id = $row["source_id"];
      $sources{$source_id}=$source_name;
   }   
   return $sources;     
}
/*************************************************************** 
   This function displays the strains count 
    for a given source id
****************************************************************/ 
function getStrainsCount($source_id){
    global $con;$straincount=0;global $SOURCEDATA;
    if($source_id<=0)return $straincount;
   if(!$con)$con=getConnection();
   $query="select count(distinct strain_id) as straincount from snp_strain_by_source ";
   $query.=" where source_id=$source_id";   
   $result = mysql_query($query,$con);
   $numRows = mysql_numrows($result);
   $row=mysql_fetch_assoc($result);
   $straincount= $row["straincount"];
   $SOURCEDATA{$source_id}{"straincount"}=$straincount; 
   //save this into memory for later use 
   return $straincount;     
}
/*************************************************************** 
   This function returns the strains count 
    for each source id
****************************************************************/ 
function loadStrainsCtBySource(){
   global $con;$straincount=0;global $SOURCEDATA;
   if(!$con)$con=getConnection();
   $query="select source_id,count(distinct strain_id) as straincount from snp_strain_by_source ";
   $query.=" group by source_id";   
   $result = mysql_query($query,$con);
   $numRows = mysql_numrows($result);
   while($row=mysql_fetch_assoc($result)){
      $straincount = $row["straincount"];
      $source_id = $row["source_id"];
      $STRAINSBYSOURCE{$source_id}=$straincount;
   }   
}

/*************************************************************** 
   This function displays the SNPs count 
    for a given source id
****************************************************************/ 
function getSNPsCount($source_id){
   global $con;$SNPcount=0;global $SOURCEDATA;
   if($source_id<=0) return $SNPcount;
   if(!$con)$con=getConnection();
   $query="select count(snpid)as snpcount from snp_by_source ";
   $query.=" where source_id=$source_id";
   $result = mysql_query($query,$con);
   $numRows = mysql_numrows($result);
   $row=mysql_fetch_assoc($result);
   $SNPcount= $row["snpcount"];
   return $SNPcount;    
}
/*************************************************************** 
   This function returns the SNPs count 
    for each source id
****************************************************************/ 
function loadSNPsCtBySource(){
   global $con;$SNPcount=0;global $SNPSBYSOURCE;
   if(!$con)$con=getConnection();
   $query="select source_id,count(snpid)as snpcount from snp_by_source ";
   $query.=" group by source_id";
   $result = mysql_query($query,$con);
   $numRows = mysql_numrows($result);
   while($row=mysql_fetch_assoc($result)){
      $SNPcount = $row["snpcount"];
      $source_id = $row["source_id"];
      $SNPSBYSOURCE{$source_id}=$SNPcount;
   }   
   
}
/*************************************************************** 
   This function returns the bp allele complement
    for the provided bp allele
****************************************************************/ 
function complementGeno($allele){
  $newallele=strtoupper($allele);
  if((strtoupper($allele)=="A"))$newallele="T";
  else if((strtoupper($allele)=="T"))$newallele="A";
  else if((strtoupper($allele)=="C"))$newallele="G";
  else if((strtoupper($allele)=="G"))$newallele="C";
  return $newallele;
  
}
//////////////////////


?>
