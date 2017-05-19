/**********************************************************************************
 * JavaScript and Ajax utils for ESTs Libraries project
 *
 * Author:  Lucie Hutchins, Scientific Software Engineer
 * Company: The Jackson Laboratory
 * Group :  Computational Sciences and System Biology
 * Creation Date :   August 2012
 *
 ********************************************************************************/
var organism;var organ; var tissue;
/*************************************************************
   This function hides or display a given column table 
   given a column, a table and the toggle option
   (0 for the X click option,1 the check box option)
**************************************************************/
function toggleColumn(columnno,tableid,type)
{
   // var this_link=document.getElementById("hq");
  var table= document.getElementById(tableid);
  var allRows = table.rows; var cellNo = -1; var j=0;
  var cells = allRows[0].cells;         
      // get the table columns header row (an array of)
      // and get the index of the column to toggle
  for(var j = 0; j< cells.length; j++){
      if(columnno.value== cells.item(j).id){ cellNo = j;}
  }
  alert(cellNo+" selected");
  for(var i=0; i<allRows.length; i++){ // hide all the row cells of this column
          if(allRows[i].cells.length > 1) {
             if(allRows[i].cells.item(cellNo).style.display=="none"){
                allRows[i].cells.item(cellNo).width="auto";
                allRows[i].cells.item(cellNo).style.display="table-cell";
	      }
	      else allRows[i].cells.item(cellNo).style.display="none";
          }
    }
}
/*****************************************
 Validates the form input data
*******************************************/
function validate(){
  if(document.getElementById("q").value==""){
     //alert("called");
     document.getElementById("wrn").innerHTML="<font color='red'>You must enter a search term</font>";
     document.getElementById("wrn").style.display="block";
  }
}
/**** for a given select list, returns the selected item's value ***/
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
 /***  clears the content of a given list ***/
 function clearSelect(thelist){
    var selLength = thelist.length;
    for(i=1; i<selLength; ++i)	{thelist.remove(1); }
 } 

/*** fill the appropriate list with the returned data from server ***/
 function fillSelect(selectListItems)
    { 
        var organismlist= selectListItems.split(",");   //list of all the organs for a given organism
	var selLength = organismlist.length;var i;
        clearSelect(document.getElementById("organism"));    // clear the previous organs list
        if(organismlist[0] !=""){
	   for(i=0; i<selLength; ++i){
               var organslist=organismlist[i].split(":");
               var org_name=organslist[1]+"("+organslist[2]+")"; 
	       addOption(document.getElementById("organism"), org_name,organslist[2]);
	    }
	 }
 }
 // add an option with 'theText' and 'theValue' to select list 'theSel'
 function addOption(theSel, theText, theValue)
    {
	var newOpt = new Option(theText, theValue);
	var selLength = theSel.length;
	theSel.options[selLength] = newOpt;
    }
function update(){
   populateOrganismList();
   updateAnnotList();
}
/***** For the selected query type , populate:
 the organism list  ****/
function populateOrganismList()
{   
   var toolid= getSelected(document.getElementById("anotselect"));  // Find the selected organism
   var url; var xmlHttp;xmlHttp= init_XmlHttpObject();
   var curDateTime = new Date(); var curHour = curDateTime.getHours();
   var curMin = curDateTime.getMinutes(); var curSec = curDateTime.getSeconds();
   var date =curHour+":"+curMin+":"+curSec;
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     url="getOrganismList.php?tool="+toolid+"&t==0&sid="+Math.random()+"&data="+date;
    // alert(url);
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if (xmlHttp.status == 200){
                           var  selectListItems;
                           var type=1;
                           // read response text into an array
                           selectListItems = xmlHttp.responseText;
                           fillSelect(selectListItems);
                       }
                       else {
                          alert('There was a problem with the request.');
                       }
                     }
          }
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime();
     var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);
     delete curDateTime;

}
function updateToolDoc(){
 var toolid=  getSelected(document.getElementById("anotselect"));
  var url; var xmlHttp;
   xmlHttp= init_XmlHttpObject();
   var curDateTime = new Date(); var curHour = curDateTime.getHours();
   var curMin = curDateTime.getMinutes(); var curSec = curDateTime.getSeconds();
   var date =curHour+":"+curMin+":"+curSec;
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     url="getOrganismList.php?tool="+toolid+"&t=1&sid="+Math.random()+"&data="+date;
     xmlHttp.onreadystatechange= function(){
        if (xmlHttp.readyState==4)
            {
                        if (xmlHttp.status == 200){
                           document.getElementById("filter").innerHTML = xmlHttp.responseText;
                       }
                       else {
                          alert('There was a problem with the request.');
                       }
                     }
          }
  xmlHttp.open("GET",url,true);
  xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
  var t = curDateTime.getTime();
  var thetimestring=t.toString(10);
  xmlHttp.send(thetimestring);
  delete curDateTime;
}
function getFile(){
var file = document.getElementById('file').files[0];
   if(file){
       blobURLref = window.URL.createObjectURL(file);
       myimg.src = blobURLref;

   }
}
/************************************************************
 For the selected organism, update the annotation listing
*************************************************************/
function updateAnnotList(){
   var orgg_id;//= getSelected(document.getElementById("orgg"));  // Find the selected organism group
   var oversion_id=getSelected(document.getElementById("organism"));  // Find the selected organism 
   var url; var xmlHttp;
   xmlHttp= init_XmlHttpObject();
   var curDateTime = new Date(); var curHour = curDateTime.getHours();
   var curMin = curDateTime.getMinutes(); var curSec = curDateTime.getSeconds();
   var date =curHour+":"+curMin+":"+curSec;
   if (xmlHttp==null)
     {
        alert ("Your browser does not support AJAX!");
        return;
     }
     url="getOrganismList.php?o="+orgg_id+"&ov="+oversion_id+"&t=1&sid="+Math.random()+"&data="+date;
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if (xmlHttp.status == 200){
                           document.getElementById("filter").innerHTML = xmlHttp.responseText;
                       }
                       else {
                          alert('There was a problem with the request.');
                       }
                     }
          }

 
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime();
     var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);
     delete curDateTime;


}
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
		  } catch (e) {alert('cant declare httpREQUEST OBJECT');}
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
  
 // select all the libraries on form
 function checkAll(){
      // boxes_count = document.organismform.libId.length;
       boxes_count = document.organismform.libId.length;
       //alert("Total boxes: "+boxes_count);
       txt = ""
       if(document.organismform.all.checked){  //select all the libraries
            for (i = 0; i < boxes_count; i++) {
                 if (!document.organismform.libId[i].checked) {
                      //txt = txt + document.f1.Liked[i].value + " "
                      document.organismform.libId[i].checked=true;
                  }
              }
        }
        else{  //unselect all the libraries
           for (i = 0; i < boxes_count; i++) {
                 if (document.organismform.libId[i].checked) {
                      //txt = txt + document.f1.Liked[i].value + " "
                      document.organismform.libId[i].checked=false;
                  }
              }
        }
 }
 
