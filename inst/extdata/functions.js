//
//
// functions.js
//
// copyright (c) 2010, Danny Arends
// last modified mrt, 2010
// first written mrt, 2010
// 
// Javascript functions for QTLviewer
//
//

function getName(evt) {
  varTextBox.firstChild.nodeValue = "Trait, Marker, LOD: " + evt.currentTarget.getAttributeNS(null,"fname");
}

function onclickName(evt) {
  createText(evt.currentTarget.getAttributeNS(null,"name"),evt.currentTarget.getAttributeNS(null,"x"),evt.currentTarget.getAttributeNS(null,"y"),getgroup(viewer.parentNode,viewer.firstChild),"12px")
}

function onCircleclickName(evt) {
  createText(evt.currentTarget.getAttributeNS(null,"name"),evt.currentTarget.getAttributeNS(null,"cx"),evt.currentTarget.getAttributeNS(null,"cy"),getgroup(viewer.parentNode,viewer.firstChild),"12px")
}

function createRect(x,y,w,h,red,green,blue,group,trait,name,lod,mqm){
  red = Math.round(red);
  green = Math.round(green);
  blue = Math.round(blue);
  if(x > 1024) return;
  if(y > 750) return;
  var newRect = document.createElementNS(svgNS,"rect");
  newRect.setAttributeNS(null,"x",x);
  newRect.setAttributeNS(null,"y",y);
  newRect.setAttributeNS(null,"width",w);
  newRect.setAttributeNS(null,"height",h);
  newRect.setAttributeNS(null,"shape-rendering","optimizeSpeed");
  newRect.setAttributeNS(null,"fill-opacity",1);
  newRect.setAttributeNS(null,"name",name+ ", " + lod + " LOD");
  newRect.setAttributeNS(null,"fname",trait+": "+name+ ", " + lod + " LOD");
  newRect.addEventListener("click",onclickName,false);
  if(mqm){
    newRect.setAttributeNS(null,"fill","rgb(255,0,0)");
  }else{
    newRect.setAttributeNS(null,"fill","rgb("+ red +","+ green+","+blue+")");
  }
  group.appendChild(newRect);
}

function createColorCircle(x,y,r,red,green,blue,group,trait,name,lod,mqm){
  red = Math.round((red+mqm*255));
  green = Math.round(green);
  blue = Math.round(blue);
  if(mqm){
    green = 0;
    blue = 0;
  }
  var newCircle = document.createElementNS(svgNS,"circle");
  newCircle.setAttributeNS(null,"cx",x);		
  newCircle.setAttributeNS(null,"cy",y);	
  newCircle.setAttributeNS(null,"r",r);
  newCircle.setAttributeNS(null,"fill","rgb("+ red +","+ green +","+ blue +")");
  newCircle.setAttributeNS(null,"stroke","rgb("+mqm*255+",0,0)");
  newCircle.setAttributeNS(null,"name",name+ ", " + lod + " LOD");
  newCircle.setAttributeNS(null,"fname",trait+": "+name+ ", " + lod + " LOD");  
  newCircle.addEventListener("mouseover",getName,false);
  newCircle.addEventListener("click",onCircleclickName,false);
  group.appendChild(newCircle);
}


function createLine(xs,ys,xe,ye,group){
  var newLine = document.createElementNS(svgNS,"line");
  newLine.setAttributeNS(null,"x1",xs);		
  newLine.setAttributeNS(null,"y1",ys);	
  newLine.setAttributeNS(null,"x2",xe);		
  newLine.setAttributeNS(null,"y2",ye);
  newLine.setAttributeNS(null,"stroke","rgb(0,0,0)");
  group.appendChild(newLine);
}

function createCircle(x,y,r,group,mqm,name,lod){
  var newCircle = document.createElementNS(svgNS,"circle");
  newCircle.setAttributeNS(null,"cx",x);		
  newCircle.setAttributeNS(null,"cy",y);	
  newCircle.setAttributeNS(null,"r",r);
  newCircle.setAttributeNS(null,"fill","rgb("+mqm*255+",0,0)");
  newCircle.setAttributeNS(null,"stroke","rgb(0,0,0)");
  newCircle.setAttributeNS(null,"name",name+ ", " + lod + " LOD");
  newCircle.addEventListener("click",onCircleclickName,false);  
  group.appendChild(newCircle);
}

function createText(mytext,x,y,group,fontsize){
  if(x > 1024) return;
  if(y > 750) return;
  var newText = document.createElementNS(svgNS,"text");
  newText.setAttributeNS(null,"x",x);		
  newText.setAttributeNS(null,"y",y);	
  newText.setAttributeNS(null,"font-size",fontsize);
  newText.setAttributeNS(null,"text-anchor","right");
  newText.setAttributeNS(null,"fill-opacity",1);		
  newText.setAttributeNS(null,"fill","rgb(1,1,1)");
  var textNode = document.createTextNode(mytext);
  newText.appendChild(textNode);
  group.appendChild(newText);
}

function deleteviewgroup(parent,child){
  while (child != null) {
    var next = child.nextSibling;

    if (child.nodeName == "group"){
      child.parentNode.removeChild(child);
    }
    child=next;
  }
}

function getgroup(paren,child){
  while (child != null) {
    var next = child.nextSibling;
    if (child.nodeName == "group"){
      return child;
    }
    child=next;
  }
}

function updateplot(){
  deleteviewgroup(viewer.parentNode,viewer.firstChild);
  if(plottype==1){
    plotInfoBox.firstChild.nodeValue = "Heatmap of "+ traits;
    drawheatmap();
  }
  if(plottype==2){
    if(locdata==null){
      plotInfoBox.firstChild.nodeValue = "Heatmap (No locations) of "+ traits;
      drawheatmap();
    }else{
      plotInfoBox.firstChild.nodeValue = "CisTrans of "+ traits;
      drawcistrans();
    }
  }
  if(plottype==3){
    plotInfoBox.firstChild.nodeValue = "Circle of 1 ";
    drawcircle();
  }
  if(plottype==4){
    plotInfoBox.firstChild.nodeValue = "Profile of 1 ";
    drawprofile();
  }
  plotInfoBox.firstChild.nodeValue += " traits examined at " + markers + " markers from " + chromosomes + " chromosomes";
}

function drawprofile(){
  if(ytrans >=traits) ytrans=traits-1;
  var ny,nx,chr,cm,mqm,ncm;
  var lodscore;
  var maxitems = absmaxitems;
  var group = document.createElementNS(svgNS,"group")
  createText(qtldata[ytrans][0],400,140,group,40);    
  var pcm=0;
  var pchr=0;
  createLine(95,140,95,705,group);
  for(var x=0;x<=maxQTL;x=x+5){
    createText(x,70,705 - (scaleG*x),group,12);    
  }
  createText("LOD score",70,140,group,12);
  for(var x=0;x+xtrans < markers;x++){
    nx = x+xtrans;
    mqm=0;
    if(modeldata!=null){
      mqm = modeldata[ytrans][nx];
    }
    lod = qtldata[ytrans][nx+1];
    chr = (mapdata[nx][1]-1);
    cm  = mapdata[nx][2];
    if(plotby==1){
      createColorCircle(100+(15*chr)+(3*x*zoomlevel),700 - (scaleG*lod),3,scaleR*lod,255-scaleR*lod,255-scaleB*lod,group,qtldata[ytrans][0],mapdata[nx][0],lod,mqm);
    }else{
      ncm  = cm + chrL[(chr)]+15*zoomlevel*(chr);
      if(x==0) pcm = ncm;
      createColorCircle(100+((ncm-pcm)*zoomlevel),700 - (scaleG*lod),3,scaleR*lod,255-scaleR*lod,255-scaleB*lod,group,qtldata[ytrans][0],mapdata[nx][0],lod,mqm);
    }
  }
  if(plotby==1){  
    createLine(95,705,100+(15*chr)+(3*x*zoomlevel),705,group);  
  }else{
    createLine(95,705,105+((ncm-pcm)*zoomlevel),705,group);  
  }
  viewer.appendChild(group);
}

function drawcistrans(){
  var ny,nx,chr,cm,ncm;
  var maxitems = absmaxitems;
  var group = document.createElementNS(svgNS,"group");
  var pcm=0;
  var mqm=0;
  for(var y=0;y<maxitemsy && (y+ytrans)<traits;y++){
    ny = y+ytrans;
    var tchr = locdata[ny][1];
    var tcm = locdata[ny][2] + chrL[(tchr-1)];    
    createText(qtldata[ny][0],0,725-((25*tchr) + (tcm)),group,"10px");    
    for(var x=0;x<maxitems && (x+xtrans) < markers;x++){
      nx = x+xtrans;
      if(modeldata!=null){
        mqm = modeldata[ny][nx];
      }
      lod = qtldata[ny][nx+1];
      chr = mapdata[nx][1]-1;
      cm  = mapdata[nx][2];
      ncm  = cm + chrL[(chr)]+15*zoomlevel*(chr);;
      if(x==0) pcm = ncm;
      if(plotby==1){      
        createRect(250+(15*chr*zoomlevel)+(3*x*zoomlevel),725-((25*tchr) + (tcm)),(2*zoomlevel)+2,10,255-scaleR*lod,255-scaleG*lod,225-scaleB*lod,group,qtldata[ny][0],mapdata[nx][0],lod,mqm);        
      }else{
        createRect(250+(2*(ncm-pcm)*zoomlevel),725-((25*tchr) + (tcm)),(2*zoomlevel)+2,10,255-scaleR*lod,255-scaleG*lod,225-scaleB*lod,group,qtldata[ny][0],mapdata[nx][0],lod,mqm);        
      }
    }
  }
  viewer.appendChild(group);
}

function drawheatmap(){
  var ny,nx,cm,chr,ncm;
  var maxitems = absmaxitems;
  var group = document.createElementNS(svgNS,"group");
  var pcm=0;
  var mqm=0;
  for(var y=0;y<maxitemsy && (y+ytrans)<traits;y++){
    ny = y+ytrans;
    createText(qtldata[ny][0],0,140+(15.15*y),group,"10px");  
    for(var x=0;x<maxitems && (x+xtrans) < markers;x++){
      nx = x+xtrans;
      if(modeldata!=null){
        mqm = modeldata[ny][nx];
      }      
      lod = qtldata[ny][nx+1];
      chr = (mapdata[nx][1]-1);
      cm  = mapdata[nx][2];      
      if(plotby==1){
        createRect(250+(15*chr*zoomlevel)+(3*x*zoomlevel),135 + (15*y),(2*zoomlevel)+2,10,255-scaleR*lod,255-scaleG*lod,225-scaleB*lod,group,qtldata[ny][0],mapdata[nx][0],lod,mqm);
      }else{
        ncm  = cm + chrL[(chr)]+15*zoomlevel*(chr);;
        if(x==0) pcm = ncm;  
        createRect(250+(2*(ncm-pcm)*zoomlevel),135 + (15*y),(2*zoomlevel)+2,10,255-scaleR*lod,255-scaleG*lod,225-scaleB*lod,group,qtldata[ny][0],mapdata[nx][0],lod,mqm);        
      }
    }
  }
  viewer.appendChild(group);
}


function drawcircle(){
  if(ytrans >=traits) ytrans=traits-1;
  var x=0;
  var y=0;
  var Px,Py;
  maxitems = absmaxitems;
  var group = document.createElementNS(svgNS,"group");
  createText(qtldata[ytrans][0],350,140,group,40);    
  var pcm=0;
  var maxchr = (mapdata[markers-1][1]);
  var mlength = chrL[maxchr];
  for(x=0;x < markers;x++){
    var mqm=0;
    if(modeldata!=null){
      mqm = modeldata[ytrans][x];
    }      
    lod = qtldata[ytrans][x+1];
    var chr = (mapdata[x][1]-1);
    if(plotby==1){
      Px = Math.sin(Math.PI*2*((x+5*chr)/(markers+5*maxchr)))
      Py = Math.cos(Math.PI*2*((x+5*chr)/(markers+5*maxchr)))
    }else{
      Px = Math.sin(Math.PI*2*((chrL[chr]+mapdata[x][2]+5*chr)/(mlength+5*maxchr)))
      Py = Math.cos(Math.PI*2*((chrL[chr]+mapdata[x][2]+5*chr)/(mlength+5*maxchr)))
    }
    createCircle(450+100*Px,400 - 100*Py,mqm*5+1,group,mqm,mapdata[x][0],lod);  
    createColorCircle(450+225*Px,400 - 225*Py,2+lod,255-scaleR*lod,255-scaleG*lod,200-scaleB*lod,group,qtldata[ytrans][0],mapdata[x][0],lod,0);
  }
  viewer.appendChild(group);
}