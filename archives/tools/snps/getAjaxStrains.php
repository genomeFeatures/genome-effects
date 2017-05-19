<?php

   $newpath= ini_get('include_path');
   $newpath .=":/raid/Users/LucieH/utils";
   $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils";
   ini_set('include_path', $newpath);
   require_once("build37/global_utils.php");
   $source_id=16; //default imputed source
   $group_id=4; //default CC Founders
   if(!empty($_GET["srcid"]))$source_id=$_GET["srcid"];
   if(!empty($_GET["gid"]))$group_id=$_GET["gid"];
   //now get this source strains list
    $strains=getThisSourceStrains($source_id,$group_id);   
   echo "$strains";
?>
