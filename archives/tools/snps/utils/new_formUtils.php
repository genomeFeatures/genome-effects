<?php

/*********************************************************** 
  This file stores 
  all functions specific to the main page
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Date : November 2009

***************************************************************/ 
 $newpath= ini_get('include_path');
  $newpath .=":/xraid/Users/LucieH/utils";
  $newpath.=":/srv/www/htdocs/cgdsnpdb/utils/utils";
  ini_set('include_path', $newpath);
  require_once("build37/global.php");

  global $username; global $passwd;global $database;global $servername;
  global $SNP_HOME; global  $SNP_SEARCH_HOME;global $CGD_HOME; global $SOURCEDATA;
   $taboutput="";$htmloutput=""; $bdconfailed=0;
 

/*************************************************************** 
   This function fills the  main strains selection
****************************************************************/   
function getCheckBoxStrainList($con)
{
     $i=1; $strains=getStrainMainList($con);
     
     $data='<table cellspacing="0" cellpadding="0" id="slist"><tr>';
     global $CCFOUNDER;
     foreach($strains as $strainID=>$strainName){   // for each strain get the list of sources
              $itemid="strain$strainID";
              $imageid="strImg$strainID";
              $displaytype="block";
              $source_list=getStrainSourceList($strainID);
              if($i==1)$data.="<td valign='top'><ul>";
              $data .="<li style='z-index:0; font-size:0.80em; padding:0;'>
                           <input type='checkbox' 
                              name='st' id='s$strainID' value='$strainID' onmouseover=\"showThisItem('$itemid','block');\" ";
              $data.="        onmouseout=\"hideThisItem('$itemid');\" onclick=\"togleImage('s$strainID');\"";
              if($CCFOUNDER{$strainID}=="$strainName")$data.=" checked='checked'";
              $data .="    />$strainName
                          <fieldset id='$itemid' style='display:none;position:absolute; 
                             padding:0; margin:0; z-index:4; background:#3F2C6D;'>
                             $source_list
                           </fieldset></li>";
               if(($i%20==0)&& ($i!=1))$data .="</ul></td><td valign='top'><ul>"; //
               if(($i==count($strains)))$data.="</ul></td>";
                $i++;
          }        
       $data .="</tr></table>";
       return $data;
}
/*************************************************************** 
   This function displays the selected strains images
****************************************************************/   
function getStrainImagesList($con)
{
     $i=1; $strains=getStrainMainList($con);
     $data="<table id='Straintable'><tr><th colspan='2' class='sectionName'>Selected Strains</th></tr>";
     
     foreach($strains as $strainID=>$strainName){   // for each strain get the list of sources
       $itemid="strain$strainID";
       $imageid="strImg$strainID";
       $displaytype="block";
       if($strainID==7){
            $data.="<tr id='$imageid'><td class='inbred'>$strainName</td>";
            $data .='<td><a href="http://phenome.jax.org/pub-cgi/phenome/mpdcgi?rtn=strains/details&strainid=664">
                        <img alt="C57BL/6J" src="http://phenome.jax.org/phenome/protodocs/Jax4/strainpix/000664_t.jpg" 
                         width="100" height="50" border="0"/></a></td></tr> '; 
       }      
       elseif($strainID==2){
           $data.="<tr id='$imageid'><td class='inbred'>$strainName</td>";
           $data.=' <tr id="t000646"><td class="inbred">A/J</td>
                         <td><a href="http://phenome.jax.org/pub-cgi/phenome/mpdcgi?rtn=strains/details&strainid=646">
                       <img alt="A/J :Inbred" src="http://phenome.jax.org/phenome/protodocs/Jax4/strainpix/000646_t.jpg"
                           width="100" height="50" border="0"/></td></tr>';
       }   
       elseif($strainID==1){
           $data.="<tr id='$imageid'><td class='inbred'>$strainName</td>";
           $data.='<td><a href="http://phenome.jax.org/pub-cgi/phenome/mpdcgi?rtn=strains/details&strainid=2448">
                                <img alt="129S1/SvImJ :Inbred" src="http://phenome.jax.org/phenome/protodocs/Jax4/strainpix/002448_t.jpg"
                           width="100" height="50" border="0"/></a></td></tr>';
       }
       elseif($strainID==8){
           $data.="<tr id='$imageid'><td class='inbred'>$strainName</td>";
           $data .='<td><a href="http://phenome.jax.org/pub-cgi/phenome/mpdcgi?rtn=strains/details&strainid=928">
                       <img alt="CAST/EiJ :wild-derived inbred"
                           src="http://phenome.jax.org/phenome/protodocs/Jax4/strainpix/000928_t.jpg"
                           width="100" height="50" border="0"/></a></td></tr>';
        } 
        else{
           $data.="<tr id='$imageid'><td class='inbred'>$strainName</td><td>Description not available</td></tr>";
        }
     }  
     $data .='<tr><td></td><td id="hc"> 
                 <a name="hc" onclick="javascript:document.getElementById(\'showStrain\').style.display=\"block\";
                      javascript:document.getElementById(\'Straintable\').style.display=\"none\";"> 
                      Hide this list &#187;</a></td></tr> </table>';     
       return $data;
}


/*************************************************************** 
   This function generates the sources list 
    for a given strain
****************************************************************/  
function getStrainSourceList($strainid)
{
   $data="<ul style='margin:0;padding:0;list-style:none;width:8em;color:#eee;'>";
   $data .="<li style='text-decoration:underline;'><b>Provided By:</b></li>";
   $sources= getTheSourceList($strainid);
   foreach($sources as $source_id=>$source_name){
        $data.=" <li style='padding:0.2em;'>$source_name</li> ";
   }        
   return "$data</ul>";
}
/*************************************************************** 
   This function displays the providers (sources)
    of all the SNPs in CGDSNPDB
****************************************************************/  
function displaySourceList($strainid)
{
   $strainid=0; global $SNPSOURCE;global $SNPSBYSOURCE;
    $data='<table id="srclist"><tr>';
   $sources= getTheSourceList($strainid); $i=1;
   $itemcount=count($sources);
   foreach($sources as $source_id=>$source_name){
        if($i==1)$data.="<td><ul>";
        $straincount = getStrainsCount($source_id);
        $snpcount=  $SNPSBYSOURCE{$source_id};
        $data .="<li  style='z-index:0;font-size:0.80em; padding:0;'>";
        //$data.="<input type='checkbox' name='source' ";
        $data.="<input type='radio' name='source' ";
        //if($source_id==21)$data .=" checked='checked' ";
        $data.=" id='sc$source_id' value='$source_id' onclick='checkAll();'/>$source_name ($straincount :$snpcount)";
        $data .="</li>";
        if(($i%2==0)&& ($i!=1)){
             $data .="</ul></td><td><ul>";
        }
       if($i==$itemcount)$data.="</ul></td></tr>";
       $i++;
   }        
 return "$data</table>";
}

?>

