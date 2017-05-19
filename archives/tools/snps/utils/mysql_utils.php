<?php

  
 /*********************************************************** 
  This file stores 
  all functions specific to the mysql public access tutorial
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : Jan 2010
  ***************************************************************/ 
  $newpath= ini_get('include_path');
  $newpath .=":/raid/Users/LucieH/utils";
  ini_set('include_path', $newpath);
  require_once("build37/global.php");

/******************************************************
  This function returns the list of all the tables
   in the selected database
 ******************************************************/
 function getTablesList(){
     $query="show tables";
     global $con_db;
     if(!$con_db)$con_db=getConnection_db();
     $result = mysql_query($query,$con_db);
     $db_list=""; $i=0;
     while($row=mysql_fetch_assoc($result)){
       $table=$row["Tables_in_cgd_snpdb"];
       if($i==0)$db_list="$table";
       else $db_list.=":$table";
       ++$i;
      }
     return $db_list;
  }
  
/******************************************************
  This function returns the schema description of a given
   table in the selected database
 ******************************************************/
 function getTableDesc($table){
     $query="desc $table";
     global $con_db;
     if(!$con_db)$con_db=getConnection_db();
     $result = mysql_query($query,$con_db);
     $table_desc=""; 
     while($row=mysql_fetch_assoc($result)){
       $column_name=$row["Field"];
       $column_type=$row["Type"];
       $column_key=$row["Key"];
       $table_desc{$column_name}="$column_type:$column_key";
      }
   return $table_desc;
 }

 function displayDBschema($table){
   global $SNP_HOME;
   $result=""; $datalist="<div id='dtlist'><h3>List of tables in cgdsnpdb </h3><ol>";
   $data="<div id='dtlist'><form name='tb' method='get' action='$SNP_HOME/utils/public_mysql.php'>
          <b>Browse tables in cgdsnpdb</b> &nbsp;<select name='tlist' id='tlist' onchange=\"getTableDesc(this);\" ";
   $tlist=split(":",getTablesList());
   foreach($tlist as $new_table){
      if($table==$new_table)$data.="<option value='$new_table' selected='selected'>$new_table</option>";
      else $data.="<option value='$new_table'>$new_table</option>";
   }
   $data.="</select><input type='submit' name='bt' value='Show Table Description'/></form> </div>";
   $result{"select"}=$data;// $result{"list"}=$datalist;
   return $result;
 }
?>
