/******************************************************************************************
 * JavaScript and jQuery utils for Graber transcript database web service
 *
 * Author:  Lucie Hutchins, Scientific Software Engineer
 * Company: The Jackson Laboratory
 * Group :  Computational Sciences and System Biology
 * Creation Date :   April 2013
 *
 *****************************************************************************************/
var jsonService_base="http://demon.jax.org/transcriptdb/web-services/index.php?format=json";
var xmlService_base="http://demon.jax.org/transcriptdb/web-services/index.php?format=xml";

/************************************************************
 The click event triggers jQuery to collect user's selection
 from the form,composes the service request url, 
 then calls getJsonResult function
*************************************************************/
$(document).ready(function(){
    $("#qsubmit").click(function(){
        var organism_version=getSelected(document.getElementById("organism"));
        var query=document.getElementById("q").value;var annot_source="";
        var url=jsonService_base+"&t="+query+"&v="+organism_version+"&isoform=1";
        if(annot_source!="")url+="&a="+annot_source;
        getJsonResult(url);
                   
     });
 });
/***********************************************************************************
 This function calls the web service to extract the list of transcripts that overlap 
 a given query.The web service returns data in json format.
 The argument to the web service are:
 t=query  => the query term
 v=organism_version
 a=annotationSource -- optional
************************************************************************************/
function getJsonResult(url){
   $.getJSON(url,function(data) {
     var table="<tr><th>transcript</th><td>chrom</td><td>strand</td><td>txStart</td><td>txEnd</td><td>exonCount</td>";
         table+="<td>cdsStart</td><td>cdsEnd</td><td>exStarts</td><td>exEnds</td><td>gene</td><td>Source</td></tr>";
      $.each(data["annotations"], function(index, val) {
          chrom=val.annotation.chrom;tx_start=val.annotation.txStart;
          tx_end=val.annotation.txEnd;strand=val.annotation.strand;transcript=val.annotation.Name;
          gene=val.annotation.Name2;cds_start=val.annotation.cdsStart;
          cds_end=val.annotation.cdsEnd;annot_source=val.annotation.source;
          ex_starts=val.annotation.exonStarts;ex_end=val.annotation.exonEnds;
          exon_count= ex_starts.split(",");
          table+="<tr><td>"+transcript+"</td><td>"+chrom+"</td><td>"+strand+"</td><td>"+tx_start+"</td>";
          table+="<td>"+tx_end+"</td><td>"+exon_count.length+"</td>";
          table+="<td>"+cds_start+"</td><td>"+cds_end+"</td><td>"+ex_starts+"</td><td>"+ex_end+"</td>";
          table+="<td>"+gene+"</td><td>"+annot_source+"</td></tr>";
      });
       document.getElementById("res").innerHTML=table;document.getElementById("res").style.display="table";
    });
}
/**** for a given select list, returns the selected item's value ***/
function getSelected(selectList){ var selectedItem; 
  for(i=0;i<selectList.length;  ++i)
  { if(selectList.options[i].selected){selectedItem = selectList.options[i].value;break;}}
  return selectedItem;
}
 

