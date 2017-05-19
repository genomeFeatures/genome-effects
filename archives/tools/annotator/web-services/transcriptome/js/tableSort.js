
$(document).ready(function(){
 //http://countrycode.org/styles/images/sortdown.gif
//http://countrycode.org/styles/images/sortup.gif
//http://countrycode.org/styles/images/sort_icon.gif

 var table=$("table.res");  var temptable = document.getElementById("tempres");
 $('#exoncount, #functsnps,#txlen,#exlen,#cdslen')
  .wrapInner('<span title="sort this column"/>')
  .each(function(){ //get the index of the column to toggle
        var th = $(this),thIndex = th.index();
        th.click(function(){ //check if temptable is not empty
        for(var r=0;r<temptable.rows.length;++r){temptable.deleteRow(r);}
         sortTable(thIndex);});   
     });
});
function sortTable(cellNo){
 var table = document.getElementById('res');var columndata=new Array();var data="";
 var temptable = document.getElementById("tempres");
 for(var i=1; i<table.rows.length; i++) { //insert all the row cells of this column into an array for sorting
     if(table.rows[i].cells.length > 1){ data=table.rows[i].cells.item(cellNo).innerHTML+":"+i; columndata.push(data);}
  }columndata.sort(function(a,b){return parseInt(a.match(/\d+:/),10)-parseInt(b.match(/\d+:/),10);});
 for(var i=0; i< columndata.length; i++) { // order rows of the table into temp
     var fields=columndata[i].split(":");var index= parseInt(fields[1].match(/\d+$/),10);
     var row=temptable.insertRow(i);
     for(var r=0;r<table.rows[index].cells.length;++r){
         var cell =row.insertCell(r);cell.innerHTML=table.rows[index].cells.item(r).innerHTML;
     }
  }
 for(var i=0; i< temptable.rows.length; i++) {table.rows[i+1].innerHTML=temptable.rows[i].innerHTML;}
}

