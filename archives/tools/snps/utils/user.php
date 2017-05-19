<?php
 /**************************************************************
  This file handle database connections
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : December 2009
  Modified : July 2010
  
  ***************************************************************/ 
$username="pup";
$passwd="puppass";

$tx_database="graber_transcriptdb";
$tx_servername="demon.jax.org";

$tx_database_p="graber_transcriptdb";
$tx_servername_p="harlequin.jax.org";

$database="cgd_snpdb";
$servername="cgd-dev.jax.org";

$database_p="cgdsnpdb";
$servername_p="cgd.jax.org";


$database_p=$database;
$servername_p=$servername;

$tx_con=mysql_connect("$tx_servername", "$username", "$passwd");
 if (!$tx_con) connectionError(mysql_error());
 if (!mysql_select_db("$tx_database", $tx_con)){
     databaseError($tx_database,mysql_error());
 }
/*$tx_con_p=mysql_connect("$tx_servername_p", "$username", "$passwd");
 if (!$tx_con_p) connectionError(mysql_error());
 if (!mysql_select_db("$tx_database_p", $tx_con_p)){
     databaseError($tx_database_p,mysql_error());
 }
*/
$con = mysql_connect("$servername", "$username", "$passwd");
 if (!$con) connectionError(mysql_error());
 if (!mysql_select_db("$database", $con)){
     databaseError($database,mysql_error());
 }

$con_p = mysql_connect("$servername_p", "$username", "$passwd");
if (!$con_p) connectionError(mysql_error());
if (!mysql_select_db("$database_p", $con_p)){
     databaseError($database_p,mysql_error());
 }
 
 
function connectionError($mysqlError){
     echo "Unable to connect to DB: " . $mysqlError;
     //exit;
 }
function databaseError($database,$mysqlError){
     echo "Unable to select $database: " . $mysqlError;
    // exit;
 }
function queryError($query,$mysqlError){
     echo "Could not successfully run query ($query): " .$mysqlError;
    // exit;
 }
function noResultFound($query){
    echo "No result found for your search. Try to check if your search term is supported\n$query\n";
    //exit;
 }
function getConnection(){
   $con = mysql_connect("$servername", "$username", "$passwd");
   if(!$con)
     $con=getConnection_p();
   return $con;
 }
 function getConnection_p(){
   $con_p = mysql_connect("$servername_p", "$username", "$passwd");
   if(!$con_p)
      $con_p=getConnection();
   return $con_p;
 }
?>
