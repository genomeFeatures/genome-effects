 /*************************************************************** 
  This is the main javacript utils
  
  Author: Lucie Hutchins, Scientific Software Engineer
  Create : December 16 2009
  Modified : July 2010
***************************************************************/

 // is this browse a version of netscape prior to version 5?
 var NS4 = (navigator.appName =="Netscape" && parseInt(navigator.appVersion) < 5);
 var curDateTime = new Date(); var curHour = curDateTime.getHours();
 var curMin = curDateTime.getMinutes(); var curSec = curDateTime.getSeconds();
 var date =curHour+":"+curMin+":"+curSec;

/*********************************************************************
 type=0 --> update current page (links and results)  after clicking the table column x
 type= 1 -->update current page (links and results) after clicking on the 
            toggle table columns check box
 type= 2 --> browsing by page number, or snp_location filter
 type= 3 -->  browsing row count or different genotype filter
**********************************************************************/
function updatePageResult(type,this_link){
   if((type!=2)&&(type!=4)){
      update_pagesLinksHref(null,type);
      this_link=document.getElementById("hq");
   }
   else{ // browsing by page number/ snp location
       var field_n=this_link.name;
       link_href=document.getElementById("hq").href;
       this_link=document.getElementById("hq");
       this_link.href=get_updatedLinkHref(link_href,field_n,type);
   }
   getAjaxData(this_link);
   document.getElementById("snptable").style.display="table";
}

/*********************************************************************
 This function uses Ajax to get the SNP listing from the server
 for the current selection
 
***********************************************************************/
function getAjaxData(this_link){
   var myLinkAttributes= this_link.href.split("?");
   var thislink_attr=myLinkAttributes[1];
   var url;var xmlHttp;
   xmlHttp= init_XmlHttpObject();
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     usid="cgd";
     url="getAjaxNewQuery.php?"+thislink_attr;//+"&sid="+Math.random()+"&data="+date;
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if(xmlHttp.status == 200){
                            document.getElementById("tresult").innerHTML= xmlHttp.responseText;                         
                         }
                      }
          }
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime(); var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);  delete curDateTime;
  
}
/*********************************************************************
 This function update each link with the current selected fields
  for every toggling event
  
*********************************************************************/
function update_pagesLinksHref(field_n,type){
 var pagelinks =document.getElementsByTagName("a");
 var page_links_count = pagelinks.length;
 for (var i = 0; i<page_links_count; ++i) {  
      var link=pagelinks[i];
      var pattern=/plinks/i;
      var link_id=link.getAttribute("id");
      if(link_id!=null){
          var new_link=get_updatedLinkHref(link.href,field_n,type);
          document.getElementsByTagName("a")[i].href=new_link; 
       } 
  } 
}
/********************************************************************
   Given a link url, this function will update the link href field to
   reflect the current field selection
*********************************************************************/
function get_updatedLinkHref(link_href,field_n,type){
     var new_link="&sg=1";
     var myLinkAttributes= link_href.split("?");
     var row_ct=getSelected(document.getElementById("rcount"));
     if(document.mainSearch.sg.checked==false)new_link="&sg=0";
     var page_number=/page/gi;var query=/query_single/i;var rows=/rcount/i;
     var source_pat=/sourc/i; var snp_location=/snp_loc/i; 
     if(type!=3){
        if(field_n!=null)
           new_link+="&"+field_n;
      }
     if(myLinkAttributes.length>0){
           var url=myLinkAttributes[0]; var link_arguments=myLinkAttributes[1];
           var this_link_fields= link_arguments.split("&"); // get the fields list and for each field
                                                           // if the checkbox is selected, add it to the link
           var boxes_count = document.mainSearch.check.length;
           var bcount=0;
           for (i= 0; i<this_link_fields.length; ++i){
                if(this_link_fields[i].match(query)!=null)new_link+="&"+this_link_fields[i];
                else if(this_link_fields[i].match(rows)!=null){new_link+="&rcount="+row_ct;}
                else if(this_link_fields[i].match(source_pat)!=null)new_link+="&"+this_link_fields[i];
                else{
                    if((type<2)){
                       if((this_link_fields[i].match(page_number)!=null))new_link+="&"+this_link_fields[i];
                       else if((this_link_fields[i].match(snp_location)!=null))new_link+="&"+this_link_fields[i];
                    }
                    else{
                        if(type==2){ // page number provided only update snp_loc
                            if((this_link_fields[i].match(snp_location)!=null))new_link+="&"+this_link_fields[i];
                         }
                        else if(type==3){ //new rowcount initiate page number
                            if((this_link_fields[i].match(page_number)!=null))new_link+="&page=0";
                            else if((this_link_fields[i].match(snp_location)!=null))new_link+="&"+this_link_fields[i];
                        }
                         
                    }
                 }
           }
           // fields to keep
           // for the rest, check the selectd field
           var pat=/=/i;
           for(var j=0;j<boxes_count;++j){   
               if(document.mainSearch.check[j].checked==true){ 
                   if(document.mainSearch.check[j].id.match(pat)!=null)
                      new_link+="&"+document.mainSearch.check[j].id;
               }
           }
      }
     new_link=url+"?genotype=1&outype=1"+new_link; 
   return new_link;
}
/*************************************************************
   This function hides or display a given column table 
   given a column, a table and the toggle option
   (0 for the X click option,1 the check box option)
**************************************************************/
function toggleColumn(columnno,tableid,type)
{
    var this_link=document.getElementById("hq");
     if(type<2){
	var table= document.getElementById(tableid);
	var allRows = table.rows; var cellNo = -1; var j=0;
	var cells = allRows[0].cells;          // get the table columns header row (an array of)
	                                       // and get the index of the column to toggle
	for(var j = 0; j< cells.length; j++){
	       if(type==0){                   // toggling using the x option
		   if(columnno.name==cells.item(j).id){
		      cellNo = j;            // now uncheck the corresponding checkbox from the output columns
		      var boxes_count = document.mainSearch.check.length;
                      for (i = 0; i < boxes_count; i++) {
                         if(document.mainSearch.check[i].value==cells.item(j).id){
                            if(document.mainSearch.check[i].checked)
                                document.mainSearch.check[i].checked=false;
                         }
                       }
                    }
	       }
	       else{
	          if(columnno.value== cells.item(j).id){
		      cellNo = j;
		   }
	       }
        }
        for(var i=0; i<allRows.length; i++) { // hide all the row cells of this column
           if(allRows[i].cells.length > 1) {
                if(allRows[i].cells.item(cellNo).style.display=="none"){
                     allRows[i].cells.item(cellNo).width="auto";
                     allRows[i].cells.item(cellNo).style.display="table-cell";
	        }
	        else allRows[i].cells.item(cellNo).style.display="none";
            }
	}
     }
    updatePageResult(type,this_link);
}

/*******************************************************************
  This function updates the strains select list 
  according to which strain source radio button was selected
 ********************************************************************/
 function checkAll(){
       boxes_count = document.mainSearch.sStrains.length;            
       var i=0; var opIndex; 
       while(i<document.mainSearch.sStr.length){          //get the index of the selected radio button
            if(document.mainSearch.sStr[i].checked){
               opIndex=i;break;
            }
           ++i;
       }
      for(i = 0; i < boxes_count; i++) { // unselect all
          if(document.mainSearch.sStrains[i].checked) {
             document.mainSearch.sStrains[i].checked=false;
          }
      }     
       if(document.mainSearch.sStr[opIndex].value==-2){   //select all Strains 
          for(i = 0; i < boxes_count; i++) {
              if(!document.mainSearch.sStrains[i].checked) {
                  document.mainSearch.sStrains[i].checked=true;
               }
           }
        }
        else if(document.mainSearch.sStr[opIndex].value==-3){   //unselect all Strains
            for(i = 0; i < boxes_count; i++) {
                 if(document.mainSearch.sStrains[i].checked) {
                    document.mainSearch.sStrains[i].checked=false;
                  }
              }
        }
        else if(document.mainSearch.sStr[opIndex].value==-1){  //only select PERLEGEN strains
            for (i = 0; i < boxes_count; i++) {
                 if (document.mainSearch.sStrains[i].checked) {
                     if(document.mainSearch.sStrains[i].value>16) document.mainSearch.sStrains[i].checked=false;
                  }
                  else{
                      if(document.mainSearch.sStrains[i].value<=16)document.mainSearch.sStrains[i].checked=true;
                  }
             }
        }
        else if(document.mainSearch.sStr[opIndex].value==-6){  //only select CC Founders strains
            for (i = 0; i < boxes_count; i++) {
                 if(document.mainSearch.sStrains[i].checked) {
                    if((document.mainSearch.sStrains[i].value!=16)&&(document.mainSearch.sStrains[i].value!=1) 
                       &&(document.mainSearch.sStrains[i].value!=2)&&(document.mainSearch.sStrains[i].value!=8) 
                       &&(document.mainSearch.sStrains[i].value!=13)&&(document.mainSearch.sStrains[i].value!=213) 
                       &&(document.mainSearch.sStrains[i].value!=214)&&(document.mainSearch.sStrains[i].value!=7))
                       document.mainSearch.sStrains[i].checked=false;
                     
                  }
                  else{
                       if((document.mainSearch.sStrains[i].value==16)||(document.mainSearch.sStrains[i].value==1) 
                       ||(document.mainSearch.sStrains[i].value==2) ||(document.mainSearch.sStrains[i].value==8) 
                       ||(document.mainSearch.sStrains[i].value==13) ||(document.mainSearch.sStrains[i].value==213) 
                       ||(document.mainSearch.sStrains[i].value==214)||(document.mainSearch.sStrains[i].value==7))
                       document.mainSearch.sStrains[i].checked=true;
                  }
             }
        }
        else if(document.mainSearch.sStr[opIndex].value==-4){  //only select Imputed strains
               getSourceStrains(16);
        }
       else{                                                  //only select Musdivarray strains
              getSourceStrains(17);
        }

 }
/**********************************************************************************
   Given the list of strains, this function clears the strains select list and 
   only checks the selected strains
***********************************************************************************/
function fillStrains(strains){
    var boxes_count = document.mainSearch.sStrains.length;
    var strain_list=strains.split(":");
    var strainscount= strain_list.length;
    for (i = 0; i < boxes_count; i++)
         document.mainSearch.sStrains[i].checked=false;
   for (i = 0; i < boxes_count; i++) {
         for(j=0;j<strainscount;++j){
             if((document.mainSearch.sStrains[i].value==strain_list[j])||
                 (document.mainSearch.sStrains[i].value==1)){
                document.mainSearch.sStrains[i].checked=true;
              }
          }
     }
}
/************************************************************************************
  for a given select list, returns the selected item's value 
*************************************************************************************/
function getSelected(selectList){
  var selectedItem;
  var selLength = selectList.length;
  for(i=0;i<selLength;  ++i)
    {
       if(selectList.options[i].selected)
	 {
	   selectedItem = selectList.options[i].value;
           break;
	  }
   }
   return selectedItem;
}
/*********************************************************************************************
 This function uses the Ajax call to get the list of strains according to the selected source
 Imputed or Musdiv array, ...
*********************************************************************************************/
function getSourceStrains(source_id)
{
   var url;
   var xmlHttp;
   xmlHttp= init_XmlHttpObject();
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     usid="cgd";
     url="getAjaxStrains.php?usid="+usid+"&srcid="+source_id+"&sid="+Math.random()+"&data="+date;
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if(xmlHttp.status == 200){
                           var strains= xmlHttp.responseText;
                           fillStrains(strains);
                         }
                      }
          }
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime(); var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);  delete curDateTime;
}
/*********************************************************************************************
 Updates the table description section with the description of the selected table
 used to display the description of each database table (Ajax)
*********************************************************************************************/
function getTableDesc(selectList)
{
  table_id=getSelected(selectList);
   var url;var xmlHttp;
   xmlHttp= init_XmlHttpObject();
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     usid="cgd";
     url="getAjaxTableDesc.php?usid="+usid+"&tid="+table_id+"&sid="+Math.random()+"&data="+date;
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if(xmlHttp.status == 200){
                            document.getElementById("tdesc").innerHTML= xmlHttp.responseText;                         
                         }
                      }
          }
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime(); var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);  delete curDateTime;
}
/****************************************************************************
 init_XmlHttpObject creates an Ajax object to be used for subsequent
 Ajax calls
****************************************************************************/
// create an http ajax object
 function init_XmlHttpObject()
 {
      var httpRequest;
      if (window.ActiveXObject) 
       { // IE
	  try {
	        httpRequest = new ActiveXObject("Msxml2.XMLHTTP");
	   }catch (e) {
		 try {
		      httpRequest = new ActiveXObject("Microsoft.XMLHTTP");
		  } catch (e) {
		      alert('cant declare httpREQUEST OBJECT');
		    }
	     }
	} 
	else if (window.XMLHttpRequest) 
	  { 	// Mozilla, Safari, ...
	     httpRequest = new XMLHttpRequest();
	     if (httpRequest.overrideMimeType) {
	         httpRequest.overrideMimeType('text/xml');
	      }
	   }
         if (!httpRequest) {
		alert('Giving up :( Cannot create an XMLHTTP instance');
			return false;
	 }
	 return httpRequest;
  }
