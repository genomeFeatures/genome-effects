 
  //----------------------------------------------//
 
  //  Created by: Romulo do Nascimento Ferreira //
 
  //  Email: romulo.nf@gmail.com                    //
 
  //----------------------------------------------//
 
  //
 
  // Drag & Drop Columns
 
  //
 
  // NOTICE: This code may be use to any purpose without further
 
  // permission from the author. You may remove this notice from the
 
  // final code, however its appreciated if you keep the author name/email
 
  // Email me if there´s something needing improvement
 
   
 
  var drag = false;
 
  document.onmouseup = release;
 
  document.onmousemove = mouseCoords;
 
   
 
  function allowDragging(tableId) {
 
  dragTable = document.getElementById(tableId)
 
   
 
  dragTableHandler = dragTable.getElementsByTagName("thead")[0] ? dragTable.getElementsByTagName("thead")[0].getElementsByTagName("tr")[0].cells : dragTable.getElementsByTagName("tr") ? dragTable.getElementsByTagName("tr")[0].cells : null
 
   
 
  maxIndex = dragTableHandler.length
 
   
 
  tableRows = dragTable.getElementsByTagName("tr");
 
     
 
      for (x=0; x<dragTableHandler.length; x++) {
 
      makeDraggable(dragTableHandler[x])
 
      dragTableHandler[x].onselectstart = function() { return false; }
 
      }
 
   
 
  dragTable.onmouseup = release;
 
  }
 
   
 
  function makeDraggable(obj) {
 
      if(!obj) return;
 
      obj.onmousedown = function(ev){
 
          if (drag == true) return
 
          captureColumnIndex(this);   
 
          createDraggedColumn(this)
 
          drag = true
 
          return false;
 
      }
 
  }
 
   
 
  function release(e) {
 
      if (drag == false) return
 
      if (!e) e=window.event
 
      if (e.target) targ = e.target
 
      else if (e.srcElement) targ=e.srcElement
 
      orderTd(targ)
 
      drag = false
 
     
 
      if (document.getElementById("drag")) {
 
      corpo = document.getElementsByTagName("body")[0];
 
      remover = document.getElementById("drag");
 
      corpo.removeChild(remover)
 
      }
 
  }
 
   
 
  function captureColumnIndex(obj) {
 
  columnIndex = obj.cellIndex
 
  return columnIndex
 
  }
 
   
 
  function orderTd(obj) {
 
  newIndex = obj.cellIndex
 
   
 
  if (newIndex == null) return
 
  if (columnIndex == newIndex) return
 
   
 
      for (x=0; x<tableRows.length; x++) {
 
      tds = tableRows[x].cells
 
      var cell = tableRows[x].removeChild(tds[columnIndex])
 
          if (newIndex >= maxIndex || newIndex + 1 >= maxIndex) {
 
          tableRows[x].appendChild(cell)
 
          }
 
          else {
 
          tableRows[x].insertBefore(cell, tds[newIndex])
 
          }
 
      }
 
  }
 
   
 
  function createDraggedColumn(obj) {
 
   
 
  draggedTable = document.createElement("table");
 
  draggedTable.id = "drag"
 
   
 
  draggedThead = document.createElement("thead");
 
  draggedTheadTr = document.createElement("tr");
 
  draggedTheadTd = document.createElement("td");
 
  draggedTheadTd.innerHTML = tableRows[0].cells[columnIndex].innerHTML
 
  draggedTheadTr.appendChild(draggedTheadTd)
 
  draggedThead.appendChild(draggedTheadTr)
 
  draggedTable.appendChild(draggedThead)
 
   
 
  draggedTbody = document.createElement("tbody"); 
 
   
 
  for (x=1; x<tableRows.length; x++) {
 
      draggedTr = document.createElement("tr");
 
      draggedTd = document.createElement("td");
 
      draggedTd.innerHTML = tableRows[x].cells[columnIndex].innerHTML
 
      draggedTr.appendChild(draggedTd)
 
      draggedTbody.appendChild(draggedTr)
 
  }
 
     
 
  draggedTable.appendChild(draggedTbody)
 
  draggedTable.style.filter = "alpha(opacity=70)";
 
  draggedTable.style.opacity = "0.7"
 
  draggedTable.style.mozOpacity = "0.7"
 
   
 
  document.getElementsByTagName("body")[0].appendChild(draggedTable)
 
  draggedTable.style.top = posY
 
  draggedTable.style.left = posX
 
  }
 
   
 
  function mouseCoords(e) {
 
  if (!e) e = window.event
 
   
 
  if(e.pageY) {posX=e.pageX; posY=e.pageY;}
 
   
 
  else if (e.clientY) {posX=e.clientX + ietruebody().scrollLeft; posY=e.clientY + ietruebody().scrollTop;}
 
   
 
      if (document.getElementById("drag")) {
 
      dragTable = document.getElementById("drag");
 
      dragTable.style.top = posY + 3 + "px"
 
      dragTable.style.left = posX + 7 + "px"
 
      }
 
     
 
  return {x:posX, y:posY}
 
  }
 
   
 
  function ietruebody(){
 
  return (document.compatMode && document.compatMode!="BackCompat")? document.documentElement : document.body
 
      }




















































































































































































































