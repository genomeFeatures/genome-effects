<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
<html xmlns="http://www.w3.org/1999/xhtml">
<head>
<meta name="keywords" content="mouse genome,Human genome,Rat genome,genome data,transcript" />
<meta name="description" content="A comprehensive genome database containing transcript data for more than 50 organisms" />
<meta http-equiv="content-type" content="text/html; charset=utf-8" />
<title>Genome Features database</title>
<?php 
   require_once("tx_utils.php");
   $annotation=0;$ch="";$op="";$oganism_version=isset($_POST["org"])? $_POST["org"] : "mm9";
   $oganism_version= $oganism_version=="0" ? "mm9":$oganism_version;
   $query=isset($_POST["q"])?filterQuery($_POST["q"]) : "";
   $annot= isset($_POST["annot"])? $_POST["annot"]: array();
 
 ?>
<link href="css/style.css" rel="stylesheet" type="text/css" media="screen" />
<link href="css/proto.css" rel="stylesheet" type="text/css" media="screen" />
<script src="http://code.jquery.com/jquery-1.9.1.min.js"></script>
<script src="http://code.jquery.com/jquery-migrate-1.1.1.min.js"></script>
</head>
<body><form   name="organismform" action="index.php" method="post" OnSubmit="return checks(this)" enctype="multipart/form-data">
<script src="js/jquery.js"></script>
<table>
<tr><td  id="main_t"><?php include("new_header.html");?></td></tr><!-- end #menu -->
 <?php 
     global $oganism_version; global $annotations;global $query;
     $data='<tr><td id="content" valign="top" style="text-align:left;background:#fff;">';
     $data.=getToolDisplay($oganism_version,$annotations,$query);
     $data.="</td></tr><tr><td id='rsum'>";
     $header="<tr><th>transcript</td><td>chrom</td><td>strand</td><td>txStart</td><td>txEnd</td><td>exonCount</td>";
     $header.="<td>cdsStart</td><td>cdsEnd</td><td>exStarts</td><td>exEnds</td><td>gene</td><td>Source</td></tr>";
     if(!empty($query))
        $rows=getJsonTxListing($query,$oganism_version,$annotations);
     //$rows= getXmlTxListing($query,$oganism_version,$annotations);
     $data.="<table class='res' id='res'>$header$rows</table></td></tr>";
     echo $data;
 ?>
</table> </form>
<!-- end #page -->
<div id="footer"><div class="data">
	<p>Copyright (c) 2013 demon.jax.org. All rights reserved. Designed by</p>
   </div>
</div>
<!-- end #footer -->
</body>
</html>
