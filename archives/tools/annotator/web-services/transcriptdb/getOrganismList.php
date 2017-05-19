<?php
 require_once("tx_utils.php");

/* given an organism  group,
  this script returns organism_name, organism_version_id list of that organism group
*/
 $orgg_id=trim($_GET["o"]);
 $u_oversion=0;$terms="";$con=0;$t=$_GET["t"];$annotation=0;
 $tool=$_GET["tool"];
 $u_oversion=$_GET["ov"];
 if($t==0)
    $terms =displayOrganismList($tool); 
 else{
    if(!empty($tool))$terms=getXmlUrl($tool);
    else $terms=displayAnnotList($u_oversion,$annotation);
 }
//output the response
 echo "$terms";
 ?>
