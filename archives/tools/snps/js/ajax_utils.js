
  // move all selected options in 'theSelFrom' to selection list 'theSelTo'
    function moveOptions(theSelFrom, theSelTo)
    {
	var selLength = theSelFrom.length;
	var selectedText = new Array();
	var selectedValues = new Array();
	var selectedCount = 0;

	var i;

	// Find the selected Options in reverse order
	// and delete them from the 'from' Select.
	for(i=selLength-1; i>=0; i--)
	{
	    if(theSelFrom.options[i].selected)
	    {
		selectedText[selectedCount] = theSelFrom.options[i].text;
		selectedValues[selectedCount] = theSelFrom.options[i].value;
		deleteOption(theSelFrom, i);
		selectedCount++;
	    }
	}

	// Add the selected text/values in reverse order.
	// This will add the Options to the 'to' Select
	// in the same order as they were in the 'from' Select.
	for(i=selectedCount-1; i>=0; i--)
	{
		addOption(theSelTo, selectedText[i], selectedValues[i]);
		// select this option
		
	}
      //  form.selectedStrains.focus();
	if(NS4) history.go(0);
    }
  function getSelectedStrains(theSelBoxName){
      var selLength = theSelBoxName.length;
      var selectedCount = 0;
      var i;

	// Find the selected strains
	for(i=0; i<selLength; i--)
	{
	    if(theSelBoxName.options[i].selected)
	    {
		alert(" selected: "+theSelBoxName.options[i].text+" and "+theSelFrom.options[i].value+"--");
	    }
	}
  }
    // move the chosen reference strain (in selection list 'repList') back
    // to the list of strains at 'searchList'
    function unSetReference(searchList, repList)
    {
        if (repList.options.length > 0)
        {
	    repList.options[0].selected = 1;
	    moveOptions (repList, searchList);
        }
    }

    // move the selected strain in 'searchList' to selection list 'repList' to
    // become the reference strain for the query
    function setReference(searchList, repList)
    {
        ct = 0;
        selIndex = -1;

        for (i = 0; i < searchList.options.length; i++)
        {
          if (searchList.options[i].selected)
          {
              ct = ct + 1;
              selIndex = i;
          }
        }

        if (ct == 0)
        {
            alert ("Must select a strain in the 'Available Strains' box " +
	        "before you can make it the Reference Strain.");
        }
        else if (ct > 1)
        {
            alert ("Must only select one strain in the 'Available Strains' " +
	        "box before you can make it the Reference Strain.");
        }
        else
        {
	    // move any existing value from the 'Reference Strain' list back
	    // to the 'List of Strains' list:

	    if (repList.options.length > 0)
	    {
	        alert ("Your previous Reference Strain (" +
	            repList.options[0].text +
		    ") is being moved back to the 'Available Strains' box.");

	        repList.options[0].selected = 1;
	        moveOptions (repList, searchList);
	    }

	    // move the new value from the 'List of Strains' list to the
	    // 'Reference Strain' list:

	    moveOptions (searchList, repList);
        }
    }
/*****************************************************************************
 move the selected strains in selection list 'selection' down one space
 if 'down' is true, or up one space if 'down' is false
*****************************************************************************/
 function moveSelected (select, down)
    {
      if (select.selectedIndex != -1) {
        if (down) {
          if (select.selectedIndex != select.options.length - 1)
            var i = select.selectedIndex + 1;
          else
            return;
        }
        else {
          if (select.selectedIndex != 0)
            var i = select.selectedIndex - 1;
          else
            return;
        }
      var swapOption = new Object();
      swapOption.text = select.options[select.selectedIndex].text;
      swapOption.value = select.options[select.selectedIndex].value;
      swapOption.selected = select.options[select.selectedIndex].selected;
      swapOption.defaultSelected = select.options
      [select.selectedIndex].defaultSelected;
      var anIndex = select.selectedIndex;
      for (var property in swapOption){
          select.options[anIndex][property] = select.options[i][property];
      }
      for (var property in swapOption)
          select.options[i][property] = swapOption[property];
      }
      var elementcount=select.options.length;
      var j=0;
    }
     //select all strains in the user selection list 
 function selectAll(select)
    {
      var elementcount=select.options.length;
      var j=0;
      for(j=0;j<elementcount;++j)
         select.options[j].selected=1;  //select this option
    }

// This function appends a new row in table pointed by table_id
function addRow(table_id,columnName,columnIndex)
{
   // each row of table_id has four columns
   var theTable=document.getElementById(table_id);
   var rowCount=theTable.rows.length; 
   var tr =theTable.insertRow(rowCount);
   var td1 = tr.insertCell(0);// var td2 = tr.insertCell(1);
   //td1.width="200px";
   td1.innerHTML = '<div id="'+columnIndex+'" style="border:solid 1px red;">'+columnName+'</div>';
   
   //td2.innerHTML = '<input type="button" onClick="showColumn('+columnIndex+');" value="Show"/>';
   theTable.style.display="table";
} 
/*********************************************************************************************
 Updates the Strains select list according to the selected source
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


// This function will delete a row pointed by rowid
// in table with id tableid
function deleteRow(row_id,tableid)
{
     var rowCount=document.getElementById(tableid).rows.length; 
     if(rowCount>0)
        document.getElementById(tableid).deleteRow(row_id);
     else
        document.getElementById(tableid).style.display="none";
}
// add an option with 'theText' and 'theValue' to selection list 'theSel'
function addOption(theSel, theText, theValue)
    {
	var newOpt = new Option(theText, theValue);var selLength = theSel.length;
	theSel.options[selLength] = newOpt;theSel.options[selLength].selected=1;  //select this option
    }

// delete the option at 'theIndex' from selection list 'theSel'
function deleteOption(theSel, theIndex)
  {
	var selLength = theSel.length;
	if(selLength>0){theSel.options[theIndex] = null;}
    }

/****************************************************************
 I will have a hidden tabble with three columns: 
 column_header,column_id,column_status
*****************************************************************/

function deleteColumn(columnno)
{
        var table= document.getElementById("snptable");
	var allRows = table.rows;var cells = allRows[0].cells;var cellNo = -1;
	for(var j = 0; j< cells.length; j++){
		if(columnno.name == cells.item(j).id)cellNo = j;
	}
       for(var i=0; i<allRows.length; i++) {
           if(allRows[i].cells.length > 1) {allRows[i].deleteCell(cellNo);}
	}

}
function fillStrains(strains){
    var boxes_count = document.mainSearch.sStrains.length;
    var strain_list=strains.split(":");
    var strainscount= strain_list.length;
    for (i = 0; i < boxes_count; i++) {
         document.mainSearch.sStrains[i].checked=false;
     }
   for (i = 0; i < boxes_count; i++) {
         for(j=0;j<strainscount;++j){
             if((document.mainSearch.sStrains[i].value==strain_list[j])||
                 (document.mainSearch.sStrains[i].value==1)){
                document.mainSearch.sStrains[i].checked=true;
              }
          }
     }
}
/*********************************************************************************************
 Updates the Strains select list according to the selected source
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
     url="utils/getAjaxStrains.php?usid="+usid+"&srcid="+source_id+"&sid="+Math.random()+"&data="+date;
     boxes_count = document.mainSearch.sStrains.length;
     xmlHttp.onreadystatechange= function(){
                    if (xmlHttp.readyState==4)
                     {
                        if(xmlHttp.status == 200){
                           var strains= xmlHttp.responseText;
                           var strain_list=strains.split(":");
                           var strainscount= strain_list.length;
                           for (i = 0; i < boxes_count; i++) {
                              if (!document.mainSearch.sStrains[i].checked) {
                                  var strain_exist=0;
                                  for(j=0;j<strainscount;++j){
                                     if(document.mainSearch.sStrains[i].value==strain_list[j]){
                                        strain_exist=1;
                                        break;
                                     }
                                   }
                                  // now unselect this strain if not in the selected source list
                                  if(strain_exist==0) document.mainSearch.sStrains[i].checked=false;
                               }
                              else{
                                  var strain_exist=0;
                                  for(j=0;j<strainscount;++j){
                                     if(document.mainSearch.sStrains[i].value==strain_list[j]){
                                        strain_exist=1;
                                        break;
                                     }
                                   }
                                   if(strain_exist==1)document.mainSearch.sStrains[i].checked=true;
                               }
                            }
                         }
                      }
          }
     xmlHttp.open("GET",url,true);
     xmlHttp.setRequestHeader("If-Modified-Since","Sat, 1 Jan 2000 00:00:00 GMT");
     var t = curDateTime.getTime(); var thetimestring=t.toString(10);
     xmlHttp.send(thetimestring);  delete curDateTime;
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
