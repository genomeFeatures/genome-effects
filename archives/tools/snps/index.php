<?php
    session_start();  
 error_reporting(0); //disable error messages 
 $newpath= ini_get('include_path');
 $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils"; 
 $newpath .=":/raid/Users/LucieH/utils";

 ini_set('include_path', $newpath);
 require_once("build37/new_formUtils.php");

  global $SNPSBYSOURCE; global $STRAINSBYSOURCE; global $SNP_SEARCH_HOME;
  if(empty($STRAINSBYSOURCE))loadStrainsCtBySource();
  global $ASSEMBLY_BUILD;global $DBVERSION;
  $title="CGD SNP Database $DBVERSION ($ASSEMBLY_BUILD) -- The Jackson Lab";
 ?> 
 <!DOCTYPE html
 PUBLIC "-//W3C//DTD XHTML 1.0 Transitional//EN"
     "http://www.w3.org/TR/xhtml1/DTD/xhtml1-transitional.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<title><?php print"$title";?></title>
<!-- jax header -->
 <meta http-equiv="Content-Language" content="en-us"/>
  <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1"/>
  <meta name="Classification" content="research"/>
  <meta name="verify-v1" content="Asn8D9XHy175bBjQrJtq2A5Yx34ANQJaiVUpdHkxmEI=" />
  <meta name="Description" content="CGDSNP is a high quality Mouse SNP database with more 
           than 66 Millions SNPs from more than 100 strains of laboratory mice."/>
  <meta name="Keywords" content="SNPs,Mouse SNPs, SNP database,SNP,Mouse, SNP dataset,SNP"/>
 <?php global $SNP_HOME; global $CGD_HOME;
  $links="<link rel=\"stylesheet\" type=\"text/css\" href=\"$SNP_HOME/css/cgd.css\" />";
  $links.="<link rel=\"shortcut icon\" href=\"$CGD_HOME/favicon.ico\" />";
  $links.="<script type=\"text/javascript\" src=\"$SNP_HOME/js/cgdjavascript.js\"></script>";
  echo $links;
 ?>
<style>
 .adv{width:98%; margin:auto;padding:0.4em;background:#eeefff;font-size:0.85em;}
 .adv .bt{background:#7e7e7e;color:#eee;}
 .adv select, .adv label{font-size:0.90em;}

</style>
<script type="text/javascript">
 <!--
   function togleImage(triggerBoxid){
      var imageid= "t"+triggerBoxid;
      if(document.getElementById(triggerBoxid).checked){
       //show the corresponding strainImage
        document.getElementById(imageid).style.display="table-row";
      } 
      else{
        document.getElementById(imageid).style.display="none";
      }
   }
   function showThisItem(itemid,displayType){
       document.getElementById(itemid).style.display=displayType;
   }
   function hideThisItem(itemid){
      document.getElementById(itemid).style.display="none";
   }
   function validates(form){
          if(form.query_single.value=="" && form.Optionfile.value=="")
            {
                alert("Please restrict your search. All the fields are blank.");
                form.main_query.focus();
                return false;
           }
     }
   //--> 
 </script>

</head>
<body>
 
  <div class="mainWindow">
    <form name='mainSearch' action='parseQuery.php?rtn=' method='get'>
       <?php
          $title="Center for Genome Dynamics Mouse SNP Database";
          displayPageHeader($title);
          displayTopNav();
      ?>
     <table id="maintable" style="width:100%;margin:auto;"><tr>
       <td id="leftp" valign="top" >
           <?php 
              displayLeftNav();
              displaySummary();
              displayAnnouncement();
           ?>
       <div id="mouse"></div></td>
      <td id="middlep" valign="top" style="width:700px;overflow:hidden;border:solid 1px #ddd;">
         <fieldset id="gannot"><legend>SNP Source &nbsp;(&nbsp;Strains: SNPs&nbsp;)</legend>
              <?php $strainid=0;$data= displaySourceList($strainid);print "$data";?>
         </fieldset>
         <fieldset class="optionnal" id="sopt" style="display:none;border:0;width:700px;margin:auto;clear:both;float:left;">
               <legend>Select Strains</legend>
               <div id="sls">
                    <input type="radio" name="sStr" id="ssa4" checked="checked" value="4" onclick="checkAll();"/>
                       <label>CC Founders</label> 
                    &nbsp;&nbsp; <input type="radio" name="sStr" id="ssa9" value="1" onclick="checkAll();"/>
                       <label>Classical Strains</label> 
                    &nbsp;&nbsp; <input type="radio" name="sStr" id="ssa" value="-2"  onclick="checkAll();"/>
                              <label>Select all</label>&nbsp;
                               <input type="radio" name="sStr" id="ssa2" value="-3"  onclick="checkAll();" />
                               <label>Unselect all</label>
                </div>
                <div id="stc" style="border:0;width:650px;margin:auto;clear:both;">
                      <?php 
                           $data=getCheckBoxStrainList($con); 
                           print "$data";
                       ?>
                </div>
           </fieldset>
           <fieldset id="requiref" style="border:0;width:700px;margin:auto;clear:both;float:left;margin-bottom:1em;">
                    <table  style="border:0;width:700px;margin:auto;clear:both;float:left;"><tr>
                     <td style='z-index:0;'><label for="q"><b>Enter Search term <a name="qhint" 
                          onmouseover="javascript:showThisItem('qhint','block');" 
                          onmouseout="hideThisItem('qhint');">(?)</a></b></label>
                         <div id="qhint" style="display:none;"><h2>You can query using the following terms:</h2>
                              <ul><li><b>Ensembl Ids : </b> gene id, transcript, protein id (ENSMUSG00000015202)</li>
                                  <li><b>MGI genes  : </b> gene id, gene name (Cnksr3) </li>
                                  <li><b>Entrez Gene  :</b> gene id(215748)</li>
                                  <li> <b>Single chromosome:bpPosition (bp or Mbp) :</b>10:3134603</li>
                                  <li> <b>Genomic range chromosome:Position1-Position2 (bp or Mbp) : </b><br/>
                                       10: 3134304 - 3227479 OR 10: 3.1 - 3.2</li>
                                  <li><b>SNP accession id :</b>NES15060003,rs37624946,or JAX00014331</li>
                              </ul>
                        </div>
                       <input type="text" size="60" name="query_single" id="mq"/>
                       <span id="showhint" style="display:none; color:red;"></span>
                       <input type="submit" name="sb" class="cgdbutton" value="Submit Query" size="30"/>
                   </td></tr></table>
            </fieldset>
            <fieldset class="optionnal" style="border:0;width:700px;margin:auto;clear:both;float:left;"><legend>Result Options</legend>
                  <table style="width:99%;margin:auto;">
                       <tr><td valign="top">
                             <div class="res-format"><span>Format :&nbsp;&nbsp;</span>
                                 <select name="outype">
                                   <option value='1' selected="selected">HTML</option>
                                   <option value='2'>CSV</option>
                                 </select>
                              
                               <br/><span>Max Rows/Page:</span>
                                  <select name="rcount">
                                   <option value="100" selected="selected">100</option>
                                   <option value="200" >200</option>
                                   <option value="300" >300</option>
                                   <option value="500" >500</option>
                                   <option value="1000">1000</option>
                                 </select>
                                <input type="hidden" name="genotype" value="1" checked="checked"/>
                             </div>
                              <div class="res-format">
                                  <span id="title">SNP Location/Annotation </span><br/>
                                  <select name="snp_loc">
                                      <option value='12' selected="selected">Select SNP Classification</option>
                                      <option value='0'>Intronic</option>
                                      <option value='10'>Exon All(Syn, NonSyn,UTR)</option>
                                      <option value='1'>Coding Synonymous</option>
                                      <option value='2'>Coding Non Synonymous</option>
                                      <!--<option value='5'>Coding Non Synonymous - Stop Gained</option>
                                      <option value='6'>Coding Non Synonymous - Stop Lost</option>
                                      <option value='7'>Coding Non Synonymous - Start Gained</option>
                                      <option value='8'>Coding Non Synonymous - Initial Met</option>-->
                                      <option value='34'>UTR</option>
                                  </select>
                              </div>
                            </td><td valign="top">
                              <div class="res-format">
                                   <span id="title">Include In Output Columns:</span><br/>
                                   <input name="vt" type="checkbox" value="1" checked="checked"/><label>Variation Type</label>
                                   <input name="gn" type="checkbox" value="1" checked="checked"/><label>Gene Name</label>
                                   <input name="fc" type="checkbox" value="1" checked="checked"/><label>Function Class</label>
                                   <br/>
                                   <input name="cp" type="checkbox" value="1" checked="checked"/><label>CpG Sites</label>
                                   <input name="ma" type="checkbox" value="1"/><label>Minor Allele%</label>
                                   <input name="mg" type="checkbox" value="1"/><label>Missing Genotype%</label>
                              </div>
                              <div class="res-format">
                                  <div style='float:left;margin:auto;'>
                                      Remove SNPs with no variation for selected strains:
                                      <input type='checkbox' name='sg' value='1' checked/>
                                  </div>
                               </div>            
                     </td></tr></table>
                 </fieldset>
        </td></tr></table>
    <?php displayFooter(); ?>
   </form></div>
  <!-- Start of StatCounter Code -->
<script type="text/javascript">
var sc_project=4184201;
var sc_invisible=1;
var sc_partition=48;
var sc_click_stat=1;
var sc_security="dcf17143";
</script>
<script type="text/javascript" src="http://www.statcounter.com/counter/counter.js"></script><noscript><div class="statcounter"><a title="click tracking" href="http://www.statcounter.com/" target="_blank"><img class="statcounter" src="http://c.statcounter.com/4184201/0/dcf17143/1/" alt="click tracking" ></a></div></noscript>
<!-- End of StatCounter Code -->
 </body></html>
