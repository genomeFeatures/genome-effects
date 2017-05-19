<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="keywords" content="mouse genome,Human genome,Rat genome,genome data,transcript" />
<meta name="description" content="A comprehensive genome database containing transcript data for more than 50 organisms" />
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Genome Features database</title>
<?php 
   require_once("tx_utils.php");
   $annotation=0;$ch="";$op="";
   $oganism_version=isset($_POST["org"])? $_POST["org"] : "mm9";
   $orgg_id=isset($_POST["orgg"])? intval($_POST["orgg"]) : 0;
   $file=isset($_POST["file"])? intval($_POST["file"]) : "";
   $query=isset($_POST["q"])? filterQuery($_POST["q"]) : "";
   $source_id=isset($_POST[""])?filterQuery($_POST["source"]) : "21";
   $query_type=isset($_POST["queryType"])? filterQuery($_POST["queryType"]) : "";
   $oganism_version= $oganism_version=="0" ? "mm9":$oganism_version;
   $annot= isset($_POST["annot"])? $_POST["annot"]: array();
 // print_r($annot);
   $annotations="";
   if(is_array($annot)){
     foreach($annot as $annotation){
         if($annotations=="")$annotations=$annotation;
         else $annotations.=",$annotation";
     }
   }
   $org_vid=getOrganismVid($oganism_version);
   $annotations=$annotations==""?getDefaultGenePrediction($orgv_id):$annotations;
   //I need to implement a javascript function that enforces the file
   // size limit requirements (5M max size, and only Type: text/plain )
   $error="";$service_base="/tmp/webservices_uploads/";$session_id="";
  //I need to make sure I implement garbage collection to clean the webservice base
  //directory and remove every file that 
  //I need to make sure the files are not overwritten unintentionally
  //meaning two users with the same input file name
  $lastMod=$_SERVER['PATH_INFO'].$_FILES["file"]["size"];
  //$eTag=$_SERVER['HTTP_IF_NONE_MATCH'];
  $eTag=base64_encode($lastMod);//
   if(!empty($_FILES["file"]["name"])){ //  5MB maximum file size
      $MAXIMUM_FILESIZE = 5 * 1024 * 1024; 
      if($_FILES["file"]["error"] > 0){ 
         $error="File upload failed.Return Code: " . $_FILES["file"]["error"] . "<br/>";}
      else{
         if($_FILES["file"]["size"]>$MAXIMUM_FILESIZE){
            $error="The file is too big. The maximum size is 5MB. Your file size is :";
            $error.=($_FILES["file"]["size"]/1024/1024)."Mb <br/>";
          }
         if($_FILES["file"]["type"]!="text/plain"){
            $error="Bad file type. Only text files are acceptable - but you provided : ".$_FILES["file"]["type"]."<br/>";
          }
        else{
            $file="/tmp/webservices_uploads/$eTag-" . $_FILES["file"]["name"];
            if(!file_exists("$file"))
            { move_uploaded_file($_FILES["file"]["tmp_name"],$file);
              chmod($file, 0644);
            }
         }
      }
   }
 echo $error;
 ?>
<link href="css/style.css" rel="stylesheet" type="text/css" media="screen" />
<link href="css/proto.css" rel="stylesheet" type="text/css" media="screen" />
 <script type="text/javascript" src="js/estlib.js"></script>
 <script type="text/javascript" src="https://ajax.googleapis.com/ajax/libs/jquery/1.9.0/jquery.min.js"></script>
</head>
<body>
<div id="main_t" width="915"><?php include("new_header.html");?></div><!-- end #menu -->
 <form   name="organismform" action="index.php" method="post" OnSubmit="return checks(this)" enctype="multipart/form-data">
<table>
<tr><td id="content" valign="top"><div class="post" style="background:#eee;"> 
      <?php 
            global $oganism_version; global $annotations;global $orgg_id;global $query;global $query_type;
            $data=getToolDisplay($oganism_version,$annotations,$query,$query_type);
            print "$data";
       ?>
  </div></td></tr>
<!-- result display starts here -->
<tr><td id="rsum">
 <div class="data" > <script type="text/javascript" src="js/tableSort.js"></script>
 <!--<div id="qsumdata"> -->
  <?php
   global $oganism_version;global $annotation;global $query;global $file; $strainIndex=array(); global $queryTypes;
   if((empty($organism_version))&& (empty($query)&&( empty($file)))){
       $data=displayMain(); print "$data";
   }
   else{ //check the query type
     $snp_query=0;$f_class="1,2,5,6,7,8";$f_name=$queryTypes{$query_type};$snp_fields="";$rows="";$total=0;
     $source_id=21; $tx_query=0; $tx_fields="";
     if($query_type=="fsnp"){$snp_query=1;}
     else if($query_type=="ssnp"){$snp_query=1;$f_class="1";$f_name="SNPs Coding Synonymous aminoacids";}
     else if($query_type=="nssnp"){$snp_query=1;$f_class="2,5,6,7,8";$f_name="SNPs Coding non Synonymous aminoacids";}
     else if($query_type=="txsnp"){$snp_query=1;$f_class="";$f_name="SNPs On Genes";}
     $data="";
     $display="<div id='snp-view'><div class='sh'>Result Format</div><ul><li>
                     <input type='radio' value='tabview' name='resultview' checked>Tabular View</li>";
     $display.="<li><input type='radio' value='graphview' name='resultview' >Graphical View</li></ul></div>";
     if($snp_query){$is_table=1;
        $snp_fields.="<th id='t1'>SNP accession(s)</th><th id='t2'>chrom:bp position</th>
                     <th id='t3'>Alleles</th>".displaySNPsStrains($source_id,$is_table,$strainIndex);
        $snp_fields.="<th>mutation type</th>";
        $snp_fields.="<th>SNP function class</th><th>Reference Codon</th><th>SNP Codon</th>
                      <th>Pos within CDS/protein</th><th>Exon rank</th><th>transcript</th>";
        $snp_fields.="<th>gene</th>";
     }
     else{
          if(($query_type=="tx")||($query_type=="exon")){$tx_query=1;$is_table=1;
               $tx_fields.="<th>chrom: Start-End: strand</th><th>transcript(s)-gene(s)</th>
                            <th>cdsStart</th><th>cdsEnd</th><th class='sort' id='exoncount'><div><span id='log'></span><span  id='tex'>exons count</span></div></th><th id='functsnps'>FunctionalSNPs</th>";
               if($query_type=="exon")$tx_fields.="<th>exonStarts</th><th>exonEnds</th>";
               else $tx_fields.="<th id='txlen'>Tx-Len(bp)</th><th id='exlen'>ExonsLen(bp)</th><th id='cdslen'>cdsLen(bp)</th>
                                  <th>Tx-%Exons</th><th>Ex-%CDS</th>";
          }
      }
     $isoform=1;$head="";
     if(!empty($file)){ 
       $handle = fopen("$file", "r"); $no_header=0;
       $headers=`head -1 $file`; $fields=split("\t",$headers);
       if(count($fields==1))$no_header=1;
       else if(count($fields>2)){
          //foreach($header as $field)$data.="<th>$field</th>";
          // if($snp_query){$data.=$snp_fields;$data="<tr>$data</tr>";}
       }
       else{$error="";}
       if($handle) { $query="File upload - ".$_FILES["file"]["name"];
          if($snp_query){$data.=$snp_fields;$data="<tr>$data</tr>";$blocks="";
              while (($line = fgets($handle)) !== false) {
                   $line=trim($line);
                   $line= filterQuery($line); 
                   $blocks.= getSNPsListing($line,$f_class,$line,$total,$source_id,$strainIndex);
                }
             $rows=$blocks;
          }
          else if($tx_query){
               $data.=$tx_fields;$data="<tr>$data</tr>";$blocks="";
               while (($line = fgets($handle)) !== false) {
                   $line=trim($line);$line= filterQuery($line);
                   $blocks.= getTxListing($line,$oganism_version,$isoform,$total,$head,$query_type,$annotations);
                }
             $rows=$blocks;
          }
         fclose($handle);
        }
     }
     else{ // a query using the term search
        $line="";
        if($snp_query){$data.=$snp_fields;$data="<tr>$data</tr>";
           $rows= getSNPsListing($query,$f_class,$line,$total,$source_id,$strainIndex);
        }
        else{
            if($tx_query){
               $data.=$tx_fields;$data="<tr>$data</tr>";
               $rows= getTxListing($query,$oganism_version,$isoform,$total,$head,$query_type,$annotations);
            }
         }
      }
      $feature_type="<div id='agg'><div class='sh'>Genome Features Display Type</div>";
      $feature_type.="<ul><li><input type='radio' name='isoform' value='1' checked/>Isoform display</li>";
      $feature_type.="<li><input type='radio' name='isoform' value='0' />Aggregate display</li></ul></div>";
      $options="<div id='seq-option'><div class='sh'>Genome Sequence display</div>";
      $options.="<ul><li><input type='radio' value='noSeq' name='seqview' checked>Hide</li>";
      $options.="<li><input type='radio' value='tx' name='seqview'>Transcript</li>";
      $options.="<li><input type='radio' value='cdsd' name='seqview'>CDS DNA</li>";
      $options.="<li><input type='radio' value='ex' name='seqview'>CDS RNA</li></ul></div>";
     $summary="<div class='total'> <table id='total'><tr><td>
                <div id='tsum'>Query: $query -Organism:$oganism_version<br/>Query Type: $f_name <br/>Total match: $total</div>$display</td><td>";
     if($snp_query){$is_table=0;
          $summary.=getSNPsSources($source_id)."".displaySNPsStrains($source_id,$is_table,$strainIndex);
      }
      else if($tx_query){
          $summary.=displayAnnotList($oganism_version,$annotations);
          $summary.="<table id='opt'><tr><td valign='top'>$feature_type</td><td valign='top'>$options</td></tr></table>";
      }
     $summary.="</td></tr>";
     if($head!="")$summary.="<tr><td colspan='2'>$head</td></tr>";
     $summary.="</table></div>";//$chrosomeFeatDensity=getSNPsFeatureDensity($query);
     $data="$summary<div id='snptable'><table class='res' id='res'>$data$rows</table><table id='tempres' style='display:none;'></table><br/><br/></div>";
     print $data;
   }
  ?>
 </div>
 </td></tr></table> </form>
<!-- end #page -->

<div id="footer"><div class="data">
	<p>Copyright (c) 2010 demon.jax.org. All rights reserved. Designed by</p>
   </div>
</div>
<!-- end #footer -->
</body>
</html>
