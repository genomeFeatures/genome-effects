<?php
/**************************************************************
  This file handle database connections
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : December 2009
  Modified : July 2010
  
  ***************************************************************/ 
$username="pup";
$passwd="puppass";
$ucsc_dbhost="genome-mysql.cse.ucsc.edu";
//mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A
//$connection="mysql -h$dbhost -u$user -A";$tm = localtime;

$tx_database="graber_transcriptdb";
$tx_servername="demon.jax.org";

$tx_database_p="graber_transcriptdb";
$tx_servername_p="demon.jax.org";

$cgd_database="cgd_snpdb";
$cgd_servername="cgd-dev.jax.org";

$cgd_database_p="cgdsnpdb";
$cgd_servername_p="cgd.jax.org";

 // define that client supports the multiple statements
 define('CLIENT_MULTI_STATEMENTS',65536);
 // define that client supports multiple results
  define('CLIENT_MULTI_RESULTS',131072);
 // the values of these defines I've found in the sourcecode of MySQL5: mysql_com.h
 // Attention: these values can be changed in other distributions of mysql

$con=mysql_connect("$tx_servername", "$username", "$passwd");//,null,CLIENT_MULTI_RESULTS);
 if (!$con) connectionError(mysql_error());
 if (!mysql_select_db("$tx_database", $con)){
     databaseError($tx_database,mysql_error());
 }
$con_p=mysql_connect("$tx_servername_p", "$username", "$passwd");
 if (!$con_p) connectionError(mysql_error());
 if (!mysql_select_db("$tx_database_p", $con_p)){
     databaseError($tx_database_p,mysql_error());
 }


/*
$cgd_con=mysql_connect("$cgd_servername_p", "$username", "$passwd");//,null,CLIENT_MULTI_RESULTS);
 if (!$cgd_con) connectionError(mysql_error());
 if (!mysql_select_db("$cgd_database", $cgd_con)){
     databaseError($cgd_database,mysql_error());
 }
$cgd_con_p=mysql_connect("$cgd_servername_p", "$username", "$passwd");
 if (!$cgd_con_p) connectionError(mysql_error());
 if (!mysql_select_db("$cgd_database_p", $cgd_con_p)){
     databaseError($cgd_database_p,mysql_error());
 }
*/
function connectionError($mysqlError){
     echo "Unable to connect to DB: " . $mysqlError;
     //exit;
 }
function databaseError($database,$mysqlError){
     echo "Unable to select $database:==$tx_servername " . $mysqlError;
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
function getTxConnection(){
   $con = mysql_connect("$tx_servername", "$username", "$passwd");
   if(!$con)
     $con=getTxConnection_p();
   return $con;
 }
 function getTxConnection_p(){
   $con_p = mysql_connect("$tx_servername_p", "$username", "$passwd");
   if(!$con_p)
      $con_p=getTxConnection();
   return $con_p;
 }
function getCgdConnection(){
   $con = mysql_connect("$cgd_servername", "$username", "$passwd");
   if(!$con)
     $con=getCgdConnection_p();
   return $con;
 }
 function getCgdConnection_p(){
   $con = mysql_connect("$cgd_servername_p", "$username", "$passwd");
   if(!$con)
      $con=getCgdConnection();
   return $con;
 }
?>
 
