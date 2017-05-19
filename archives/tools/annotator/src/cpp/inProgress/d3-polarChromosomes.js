//polarChromosomes.js
//Justin Hendrick
//Justin.Hendrick@jax.org
//JustinJHendrick@gmail.com
//started on 6-6-13
//graph R/CAPE results with D3.js

/*------------------------------TODO--------------------------------
    fish-eye
        low priority
        keyboard toggle or hold
    more modular
        splitting into multiple files is more hassle than it's worth
        and slower to load multiple files anyway
    selection control with keyboard
    effectPlot
        bug: error bar exceeds box
        hover tells number of points
        low priority: option for other type of effPlot
    upload data or R/clickme
        go local instead of server
    zoom & drag
        touch screen
            double tap distance threshold
    threshold slider
------------------------------------------------------------------*/
"use strict";
//read in data
d3.json("capeOut_min.json", main);

//main is a callback function. It is run after loading is done
function main(capeOut) {

    ///////////////////////////////FUNCTIONS///////////////////////////////////////////////
    function greenOrRed(d, i) {//green is positive, red is negative. This is a global color scheme
        if (d && d.eff > 0) {
            return "green"
        } else {
            return "red"
        }
    }

    function sumTo(arr, to) { //does not include (to)th element
        var s = 0;            //sum(arr) is same as sumTo(arr, arr.length)
        for(var i = 0; i < to; i++) {
            s += arr[i];
        }
        return s;
    }

    function markerPos(d) { //given a marker, returns location on circle in polar coordinates [r, a]
        var angle = theta(sumTo(chromoLens, d.chromosome) + d.cM * .01 * (chromoLens[d.chromosome] + chromoPad * totalBP / (2 * Math.PI)));
        if (isNaN(angle)) {//covariate or phenotype
            return [NaN, NaN];
        }
        return [innerRad, angle];
    }

    function fadeOthers(dat, ind) { //select all paths and polygons, filter out selected, set to low or 0 opacity depending on type.
        if(!hold) {
            ops = [];
            var connChromos = [ind];
            var others = svg.selectAll("path, polygon")
                .filter(function(d, i) {
                    var thisID = d3.select(this).attr("id");
                    var thisClass = d3.select(this).attr("class");
                    //if(!d) {
                    //    return true;
                    // }
                    if(d.target || d.source) {
                        //this is a link. Either Geno or pheno
                        if (ind == d.source.chromosome) {
                            connChromos.push(d.target.chromosome);
                            return false;
                        }
                        if (ind == d.target.chromosome) {
                            connChromos.push(d.source.chromosome);
                            return false;
                        }
                    }
                    if (thisClass != "arrows" && thisClass != "arrowHeads" && connChromos.indexOf(parseInt(thisID.slice(thisID.length - 2, thisID.length))) != -1) {
                        return false;
                    }
                    return true
                })
                .attr("opacity", function(d, i) {
                    var selThis = d3.select(this);
                    if (selThis.attr("class") == "arrows" || selThis.attr("class") == "arrowHeads") {
                        selThis.attr("visibility", "hidden");
                        return selThis.attr("opacity");
                    } else {
                        ops[i] = selThis.attr("opacity");//store the opacity in the ops array
                        return .1 * ops[i];
                    }
                })
            ;
        }
    }

    function sharpOthers(dat, ind) { //select all paths and polygons, filter out selected, bring back to original opacity
        if(!hold) {
            var connChromos = [ind];
            var others = svg.selectAll("path, polygon")
                .filter(function(d, i) {
                    var thisID = d3.select(this).attr("id");
                    var thisClass = d3.select(this).attr("class");
                    //if(!d) {
                    //    return true;
                    // }
                    if(d.target || d.source) {
                        //this is a link. Either Geno or pheno
                        if (ind == d.source.chromosome) {
                            connChromos.push(d.target.chromosome);
                            return false;
                        }
                        if (ind == d.target.chromosome) {
                            connChromos.push(d.source.chromosome);
                            return false;
                        }
                    }
                    if (thisClass != "arrows" && thisClass != "arrowHeads" && connChromos.indexOf(parseInt(thisID.slice(thisID.length - 2, thisID.length))) != -1) {
                        return false;
                    }
                    return true
                })
                .attr("opacity", function(d, i) {
                    if(ops[i]) {
                        return ops[i];
                    } else {
                        return d3.select(this).attr("opacity");
                    }
                })
                .attr("visibility", "visible")
            ;
        }
    }

    function toggleHold(d, i) { //turn chromosomal hold off and on
        var mouse = d3.mouse(document.getElementById("irisSvg"));
        function compare(a, b) { //sort by centimorgan position (increasing)
            if (a.source.cM < b.source.cM)
                return -1;
            if (a.source.cM > b.source.cM)
                return 1;
            return 0;
        }
        if(Math.pow(chromoDown[0] - mouse[0], 2) + Math.pow(chromoDown[1] - mouse[1], 2) < 25) { //if drag dist < 5, it counts as a chromo click
            if (hold && holdOn == i) { //only cancel hold when the held is clicked
                holdOn = -1; //let go of current selection
                hold = false; //no longer holding anything
                d3.select(this).attr("fill", "white"); //turn it back to white
                d3.selectAll("#cartBarSvg")
                    .remove();
                if(touch) {
                    sharpOthers(d, i)
                }
            } else if(!hold) { //set new hold
                d3.select(this).attr("fill", "#9090F0"); //turn it blue-ish
                holdOn = i; //hold on to this one
                if(touch) {
                    fadeOthers(d, i);
                }
                hold = true;
                var phenoDatOne = []; //pheno data only for the selected chromosome
                for(var j = 0; j < phenoDat.length; j++) {
                    if(phenoDat[j].source.chromosome == i) {
                        phenoDatOne.push(phenoDat[j]);
                    }
                }
                phenoDatOne.sort(compare); //call sort with my compare function
                cartesianBars(phenoDatOne); //spawn supplementary graphs
            } else {
                hold = false; //no longer holding anything
                d3.select("#mouseRegion0" + holdOn).attr("fill", "white"); //turn it back to white
                sharpOthers(d, holdOn);
                fadeOthers(d, i);
                d3.selectAll("#cartBarSvg")
                    .remove();
                d3.select(this).attr("fill", "#9090F0"); //turn it blue-ish
                holdOn = i; //hold on to this one
                hold = true;
                var phenoDatOne = []; //pheno data only for the selected chromosome
                for(var j = 0; j < phenoDat.length; j++) {
                    if(phenoDat[j].source.chromosome == i) {
                        phenoDatOne.push(phenoDat[j]);
                    }
                }
                phenoDatOne.sort(compare); //call sort with my compare function
                cartesianBars(phenoDatOne); //spawn supplementary graphs
            }
        }
    }

    function cartesianBars(data) { //draw effect bar charts
        var widthBarSvg = svgSize * .33;
        var heightBarSvg = svgSize;
        var barSvg = d3.select("#cartBar").append("svg")
            .attr("id", "cartBarSvg")
            .attr("width", widthBarSvg)
            .attr("height", heightBarSvg)
        ;
        var maxBarHeight = heightBarSvg / phenos.length - vChartPad;
        var cartBarHeight = d3.scale.linear()
            .domain([minEffBar, maxEffBar])
            .range([vChartPad, maxBarHeight])
        ;
        var reverseCartBarHeight = d3.scale.linear()
            .domain([maxBarHeight - vChartPad, 0])
            .range([minEffBar, maxEffBar])
        ;
        
        //draw stacked axes
        var cartAxis = barSvg.selectAll("rect.cartAxis") //xAxis
            .data(phenos)
          .enter().append("rect")
            .attr("class", "cartAxis")
            .attr("id", function(d, i) {return "cartAxis0" + i;})
            .attr("x", hChartPad)
            .attr("y", function(d, i) {return i * heightBarSvg / phenos.length;})
            .attr("width", widthBarSvg - hChartPad - 1)
            .attr("height", maxBarHeight)
            .attr("fill", "none")
            .attr("stroke", "black")
        ;
            
        var cartPhenoLabels = barSvg.selectAll("text.cartPhenoLabels") //phenotype labels going up on the left side
            .data(phenos)
          .enter().append("text")
            .attr("class", "cartPhenoLabels")
            .attr("id", function(d, i) {return "cartPhenoLabel0" + i;})
            .attr("x", textSize)
            .attr("y", function(d, i) {return .5 * maxBarHeight + (phenos.length - i - 1) * (maxBarHeight + vChartPad);})//xAxis plus some
            .attr("font-size", textSize)
            .attr("transform", function(d, i) {return "rotate(-90," + textSize + "," + d3.select("#cartPhenoLabel0" + i).attr("y") + ")"})
            .attr("text-anchor", "middle")
            .text(function(d, i) {return d + " Effect";})
        ;
        
        //line up marker to make comparing easier
        var markers = [];
        for(var i = 0; i < data.length; i++) {//need to know number of markers for splitting
            if(markers.indexOf(data[i].source.name) == -1) { //new marker
                markers.push(data[i].source.name);
            }
        }
        
        var splitPheno = [];
        for(var p = 0; p < phenos.length; p++) { //fill splitPheno with empty arrays
            splitPheno[p] = [];
        }
        for(var m = 0; m < markers.length; m++) {
            for(var i = 0; i < data.length; i++) { //isolate just the data for each pheno
                for(var p = 0; p < phenos.length; p++) {
                    if(phenos[p] == data[i].target.name && markers[m] == data[i].source.name) {
                        splitPheno[p][m] = data[i]
                    }
                }
            }
        }
        /*var maxM = 0;
        for(var p = 0; p < splitPheno.length; p++) {//find max number of markers
            len = splitPheno[p].length; 
            if(len > maxM) {
                maxM = len;
            }
        }*/
        
        //console.log(splitPheno);
        
        for(var p = 0; p < phenos.length; p++) { //draw ticks, bars, and labels for each pheno
            //draw y axis ticks
            for(var t = 0; t <= maxBarHeight + .01; t += (maxBarHeight - vChartPad) / 4) {//6 ticks
                barSvg
                    .append("rect")
                    .attr("class", "effTicks" + t)
                    .attr("id", function(d, i) {return "effTick" + t + ",0" + i;})
                    .attr("x", hChartPad - 3) //center it on y axis
                    .attr("width", 7) //one pixel overlaps axis
                    .attr("y", parseFloat(barSvg.select("#cartAxis0" + p).attr("y")) + t)
                    .attr("height", 1)
                ;
            
                //label y axis ticks with |effect|
                barSvg
                    .append("text")
                    .attr("class", "effTickLabels")
                    .attr("id", function(d, i) {return "effTickLabel" + t + ",0" + i;})
                    .attr("x", .40 * hChartPad) //move far enough for label to be on left
                    .attr("y", function(d, i) {
                        if(t == 0 & p == 0) {//the topmost label was being cut off. This moves it down a bit
                            return parseFloat(barSvg.select("#cartAxis0" + p).attr("y")) + t + 2 * textSize / 3;
                        } else {//centered on tick
                            return parseFloat(barSvg.select("#cartAxis0" + p).attr("y")) + t + textSize / 3;
                        }
                    })
                    .attr("font-size", textSize)
                    .text(reverseCartBarHeight(t).toFixed(2)) //2 decimal places
                ;
            }
            
            var barWidth = (widthBarSvg - hChartPad - 2) / markers.length - 1;
            var bars = barSvg.selectAll("rect.bars" + p)//add bars to axes
                .data(splitPheno[p])
              .enter().append("rect")
                .attr("class", "bars" + p)
                .attr("id", function(d, i) {return "bar" + p + "0" + i;})
                .attr("x", function(d, i) {return i * (widthBarSvg - hChartPad - 2) / markers.length + hChartPad + 2;})
                .attr("y", function(d, i) {
                    if(d) {
                        return (phenos.length - p) * maxBarHeight + (phenos.length - p - 1) * vChartPad - cartBarHeight(Math.abs(d.eff));
                    } else {
                        return (phenos.length - p) * maxBarHeight + (phenos.length - p - 1) * vChartPad;
                    }
                })
                .attr("width", barWidth)
                .attr("height", function(d, i) {
                    if(d) {
                        return cartBarHeight(Math.abs(d.eff));
                    } else {
                        return 0;
                    }   
                })
                .attr("fill", greenOrRed)
                .attr("stroke", "none")
                .attr("opacity", function(d, i) {
                    if(d) {
                        return sigToAlp(d.sig);
                    } else {
                        return 0;
                    }
                })
                .append("svg:title")
                .text(function(d, i) {
                    if(d) {
                        return "eff: " + d.eff + "\nsig: " + d.sig + "\nqval: " + d.qval;
                    } else {
                        return "";
                    }
                })
            ;
                
            var markerLabels = barSvg.selectAll(".markerLabels" + p)//label markers on each bar
                .data(splitPheno[p])
              .enter().append("text")
                .attr("class", "markerLabels" + p)
                .attr("id", function(d, i) {return "markerLabel" + p + "0" + i;})
                .attr("x", function(d, i) {
                    if (barSvg.select("#bar" + p + "0" + i).attr("height") > barWidth) {//taller than wide
                        //middle of bar in x dir
                        return parseFloat(barSvg.select("#bar" + p + "0" + i).attr("x")) + barWidth / 2 - textSize / 3;
                    } else {//wider than tall
                        return parseFloat(barSvg.select("#bar" + p + "0" + i).attr("x")) + barWidth / 2;
                    }
                })
                .attr("y", function(d, i) {
                    if (barSvg.select("#bar" + p + "0" + i).attr("height") > barWidth) {//taller than wide
                        //top of bar
                        return parseFloat(barSvg.select("#bar" + p + "0" + i).attr("y"));
                    } else {//wider than tall
                        //top of text touches top of bar
                        return parseFloat(barSvg.select("#bar" + p + "0" + i).attr("y")) + textSize;
                    }
                })
                .attr("text-anchor", function(d, i) {
                    if (barSvg.select("#bar" + p + "0" + i).attr("height") > barWidth) {//taller than wide
                        return "start";
                    } else {//wider than tall
                        return "middle";
                    }
                })
                .attr("fill", "white")
                .attr("font-size", textSize)
                .attr("transform", function(d, i) {
                    if (barSvg.select("#bar" + p + "0" + i).attr("height") > barWidth) {//taller than it is wide --> vertical text
                        var thisText = barSvg.select("#markerLabel" + p + "0" + i);
                        return "rotate(" + 90 + "," + thisText.attr("x") + "," + thisText.attr("y") + ")";
                    }
                })
                .text(function(d, i) {
                    if(d) {
                        return " " + d.source.name;
                    } else {
                        return "";
                    }
                })
            ;
        }
    }

    function toggleArrowPlots(d, i) { //control rendering of on-arrow-click plots
        var currentPos = d3.mouse(document.getElementById("irisSvg"));
        if(Math.pow(currentPos[0] - downPos[0], 2) + Math.pow(currentPos[1] - downPos[1], 2) < 25) { //if drag distance < 5, it counts as an arrow click
            if(arrowPlots && arrow == i) {//remove plots
                d3.selectAll("#effPlotSvg").remove();
                d3.selectAll("#netSvg").remove();
                arrowPlots = false;
                arrow = -1;
                d3.select(this)
                    .attr("stroke", greenOrRed)
                    .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
                ;
                d3.select("#arrowHead0" + i)
                    .attr("fill", greenOrRed)
                    .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
                ;
            } else if(!arrowPlots) {//draw new plots
                //console.log(d.source.name + ", " + d.target.name)
                arrow = i;
                effectPlots(d);
                miniNetwork(d);
                arrowPlots = true;
                d3.select(this)//turn this arrow and its head black.
                    .attr("stroke", "black")
                    .attr("opacity", 1)
                ;
                d3.select("#arrowHead0" + i)
                    .attr("fill", "black")
                    .attr("opacity", 1)
                ;
            } else {//switch to different arrow
                d3.selectAll("#effPlotSvg").remove();
                d3.selectAll("#netSvg").remove();
                d3.select("#arrow0" + arrow)
                    .attr("stroke", greenOrRed)
                    .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
                ;
                d3.select("#arrowHead0" + arrow)
                    .attr("fill", greenOrRed)
                    .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
                ;
                arrow = i;
                effectPlots(d);
                miniNetwork(d);
                d3.select(this)//turn this arrow and its head black.
                    .attr("stroke", "black")
                    .attr("opacity", 1)
                ;
                d3.select("#arrowHead0" + i)
                    .attr("fill", "black")
                    .attr("opacity", 1)
                ;
            }
        }
    }

    function effectPlots(link) { //draw effect plots
        var effColors = ["blue", "#FF6600"/*orange*/, "teal", "purple", "yellow", "brown"]; //I think 6 colors should be enough
        var nRow = link.effPlot[0].dat[0].length;//nRow and nCol of effect plot data
        var nCol = link.effPlot[0].dat.length;
        var effPlotWidth = .33 * svgSize;
        var effPlotHeight = svgSize;
        var padRight = .08 * svgSize; //padding on the right side of effect plots for labeling
        var maxBoxHeight = effPlotHeight / phenos.length - vChartPad;
        var effPlotSvg = d3.select("#effectPlot").append("svg")
            .attr("id", "effPlotSvg")
            .attr("width", effPlotWidth)
            .attr("height", effPlotHeight)
        ;
        //draw stacked axis boxes
        var axisBox = effPlotSvg.selectAll("rect.axisBox")
            .data(phenos)
          .enter().append("rect")
            .attr("class", "axisBox")
            .attr("id", function(d, i) {return "axisBox0" + i;})
            .attr("x", hChartPad)
            .attr("y", function(d, i) {return i * effPlotHeight / phenos.length;})
            .attr("width", effPlotWidth - hChartPad - padRight)
            .attr("height", maxBoxHeight)
            .attr("fill", "lightgrey")
            .attr("stroke", "black")
        ;
        
        for(var p = 0; p < phenos.length; p++) {
            var x = d3.scale.linear()//x positions based on num columns in dat
                .domain([0, link.effPlot[0].dat.length - 1])
                .range([hChartPad + .05 * effPlotWidth, .95 * effPlotWidth - padRight])
            ;
            var y = d3.scale.linear()//y scale to be added to axis position
                .domain([link.effPlot[p].ylim[0], link.effPlot[p].ylim[1]])
                 //inverted so that point data can be fed in directly
                .range([(phenos.length - p - 1) * effPlotHeight / phenos.length + maxBoxHeight - 2, (phenos.length - p - 1) * effPlotHeight / phenos.length + 2])
            ;
            
            //put scale behind first so lines show up on top
            for(var i = link.effPlot[p].ylim[0]; i <= link.effPlot[p].ylim[1] + .01; i += (link.effPlot[p].ylim[1] - link.effPlot[p].ylim[0]) / 5) {//scale lines and labels
                effPlotSvg.append("text")//scale labels
                    .attr("class", "effect levels")
                    .attr("x", hChartPad * .4)
                    .attr("y", function() {
                        if(Math.abs(i - link.effPlot[p].ylim[1]) < .01 && p == phenos.length - 1) {
                            return y(i) + 2 * textSize / 3;
                        } else {
                            return y(i) + textSize / 3;
                        }
                    })
                    .attr("font-size", textSize)
                    .text(i.toFixed(2))
                ;
                if (i > link.effPlot[p].ylim[0] && i < link.effPlot[p].ylim[1] - .01) {//scale lines
                    effPlotSvg.append("rect")
                        .attr("class", "scaleLine")
                        .attr("x", hChartPad)
                        .attr("y", y(i))
                        .attr("width", effPlotWidth - hChartPad - padRight)
                        .attr("height", 1)
                        .attr("fill", "darkgrey")
                    ;
                }
            }
            
            for(var r = 0; r < nRow; r++) {//one row is one line
                var effectLines = effPlotSvg.selectAll("#effectLines0" + r + "0" + p)//draw effect lines
                    .data([link.effPlot[p]])//data is an array with one element for this phenotype
                  .enter().append("path")
                    .attr("class", "effectLines")
                    .attr("id", function(d, i) {return "effectLines0" + p;})
                    .attr("d", function(d, i) {
                        var out = "";
                        for(var c = 0; c < nCol; c++) {//slope lines
                            if(c == 0) {
                                out += "M";
                            } else {
                                out += "L";
                            }
                            out += x(c) + "," + y(d.dat[c][r]);
                        }
                        for(var c = 0; c < nCol; c++) {//error bars
                            out += "M" + x(c) + "," + y(d.dat[c][r] + d.errors[r][c]) + "L" + x(c) + "," + y(d.dat[c][r] - d.errors[r][c]);
                            //ticks on ends of error bars
                            out += "M" + (x(c) - 5) + "," + y(d.dat[c][r] + d.errors[r][c]) + "L" + (x(c) + 5) + "," + y(d.dat[c][r] + d.errors[r][c]);
                            out += "M" + (x(c) - 5) + "," + y(d.dat[c][r] - d.errors[r][c]) + "L" + (x(c) + 5) + "," + y(d.dat[c][r] - d.errors[r][c]);
                        }
                        return out;
                    })
                    .attr("fill", "none")
                    .attr("stroke", effColors[r])
                    .attr("stroke-width", .005 * effPlotWidth)
                ;
                for(var c = 0; c < nCol; c++) {
                    var effectPoints = effPlotSvg.selectAll("#effectPoints0" + r + "0" + p) //draw point for mousover purposes
                        .data([link.effPlot[p]])//data is an array with one element for this phenotype
                      .enter().append("circle")
                        .attr("class", "effectLines")
                        .attr("id", function(d, i) {return "effectLines0" + p;})
                        .attr("cx", x(c))
                        .attr("cy", function(d, i) {return y(d.dat[c][r]);})
                        .attr("r", .01 * effPlotWidth)
                        .attr("fill", effColors[r])
                        .attr("stroke-opacity", 0)
                        .attr("stroke-width", 15)
                        .attr("pointer-events", "all")
                      .append("svg:title")
                        .text(function(d, i) {return "eff: " + d.dat[c][r] + "\nse: " + d.errors[r][c]})
                    ;
                }
                
                effPlotSvg.append("text") //xAxis labels on bottom
                    .attr("class", "xEffectLabels")
                    .attr("x", x(r))
                    .attr("y", (phenos.length - p - 1) * effPlotHeight / phenos.length + maxBoxHeight + textSize)
                    .attr("font-size", textSize)
                    .attr("text-anchor", "middle")
                    .text(r * (1 / nRow))
                ;
                
                effPlotSvg.append("rect") //marker2 legend
                    .attr("class", "m2Legend")
                    .attr("x", effPlotWidth - padRight + 3)
                    .attr("y", (phenos.length - p - 1) * effPlotHeight / phenos.length + textSize * (r + 1))
                    .attr("width", .4 * (padRight - 3))
                    .attr("height", 1)
                    .attr("fill", effColors[r])
                ;
                
                effPlotSvg.append("text") //marker2 legend label
                    .attr("class", "m2LegendLabel")
                    .attr("x", effPlotWidth - .5 * padRight)
                    .attr("y", (phenos.length - p - 1) * effPlotHeight / phenos.length + textSize * (r + 1) + textSize / 3)
                    .attr("font-size", textSize)
                    .text(r * (1 / nRow))
            }
            effPlotSvg.append("text")//marker1 name on bottom
                .attr("class", "m1Label")
                .attr("x", .5 * effPlotWidth)
                .attr("y", (phenos.length - p - 1) * effPlotHeight / phenos.length + maxBoxHeight + 2 * textSize)
                .attr("font-size", textSize)
                .attr("text-anchor", "middle")
                .text(link.source.name)
            ;
            
            effPlotSvg.append("text")//marker2 name on right, going down
                .attr("class", "m2Label")
                .attr("x", effPlotWidth - padRight + textSize / 3)
                .attr("y", (phenos.length - p - 1) * effPlotHeight / phenos.length + maxBoxHeight)
                .attr("font-size", textSize)
                .attr("text-anchor", "end")
                .attr("transform", "rotate(90," + (effPlotWidth - padRight + textSize / 3) + "," + ((phenos.length - p - 1) * effPlotHeight / phenos.length + maxBoxHeight) + ")")
                .text(link.target.name)
            ;
            
            effPlotSvg.append("text")//pheno label on left going up
                .attr("class", "effPlotPhenoLabel")
                .attr("x", textSize)
                .attr("y", .5 * maxBoxHeight + (phenos.length - p - 1) * (maxBoxHeight + vChartPad))
                .attr("text-anchor", "middle")
                .attr("font-size", textSize)
                .attr("transform", "rotate(-90," + textSize + "," + (.5 * maxBoxHeight + (phenos.length - p - 1) * (maxBoxHeight + vChartPad)) + ")")
                .text(phenos[p])
            ;
        }
    }

    function miniNetwork(link) { //draw network with 2 markers and all phenotypes
        var netWidth = .33 * svgSize;
        var netHeight = svgSize;
        var netSvg = d3.select("#net").append("svg")
            .attr("id", "netSvg")
            .attr("width", netWidth)
            .attr("height", netHeight)
        ;
        var markers = [link.source.name, link.target.name];
        
        //if doubly connected, find other link. put both in mLinks
        var mLinks = [link];
        if(link.dc) {
            for(var i = 0; i < links.length; i++) {
                if(links[i].dc && link.source.name == links[i].target.name && link.target.name == links[i].source.name) {
                    mLinks.push(links[i]);
                }
            }
        }
        
        //for source and target, if an effect on a phenotype exists, add it to pLinks
        var pLinks = [];
        for(var i = 0; i < markers.length; i++) {
            for (var j = 0; j < phenoDat.length; j++) {
                if(phenoDat[j].source.name == markers[i]) {
                    pLinks.push(phenoDat[j]);
                }
            }   
        }
        
        var mAndP = mLinks.concat(pLinks);
        //for each node, draw an ellipse
        var nodes = markers.concat(phenos);
        var phenoPos = d3.scale.linear()
            .domain([2, nodes.length])
            .range([netHeight, 0])
        ;
        
        var ellipses = netSvg.selectAll("ellipse.node")//draw ellipse
            .data(nodes)
          .enter().append("ellipse")
            .attr("class", "node")
            .attr("id", function(d, i) {return d;})
            .attr("cx", function(d, i) {
                if(i < 2) {
                    return .15 * netWidth;
                } else {
                    return .85 * netWidth;
                }
            })
            .attr("cy", function(d, i) {
                if(i < 2) {
                    return netHeight - (2 * i + 1) / 4 * netHeight;
                } else {
                    return phenoPos(i) - netHeight / (2 * phenos.length);
                }
            })
            .attr("rx", .1 * netWidth)
            .attr("ry", function(d, i) {
                if(i < 2) {
                    return netHeight / 8;
                } else {
                    return netHeight / (2 * phenos.length) - 10;
                }
            })
            .attr("fill", "none")
            .attr("stroke", "black")
        ;
        
        var netLabels = netSvg.selectAll("text.netLabel") //label the nodes
            .data(nodes)
          .enter().append("text")
            .attr("class", "netLabel")
            .attr("x", function(d, i) {return parseFloat(netSvg.select("#" + d).attr("cx")) - textSize / 3;})
            .attr("y", function(d, i) {return netSvg.select("#" + d).attr("cy");})
            .attr("font-size", textSize)
            .attr("text-anchor", "middle")
            .attr("transform", function(d, i) {return "rotate(90," + (netSvg.select("#" + d).attr("cx") - textSize / 3) + "," + netSvg.select("#" + d).attr("cy") + ")";})
            .text(function(d, i) {return d;})
        ;
        
        function getOnEll(nodeName, angle, rOffset, aOffset) {//given a marker name, angles, and offsets. Returns a point on the ellipse with the offsets applied
            var node = netSvg.select("#" + nodeName);
            return [parseFloat(node.attr("cx")) + (parseFloat(node.attr("rx")) + rOffset) * Math.cos(angle + aOffset), parseFloat(node.attr("cy")) - (parseFloat(node.attr("ry")) + rOffset) * Math.sin(angle + aOffset)];
        }
        
        function calcAngle(sourceName, targetName, dc) {
            if(phenos.indexOf(targetName) == -1) {//both are markers
                if(parseFloat(netSvg.select("#" + sourceName).attr("cy")) > parseFloat(netSvg.select("#" + targetName).attr("cy"))) {//up arrow
                    if(dc) {//doubly connected. Move both from center
                        return [Math.PI / 3, 5 * Math.PI / 3];
                    } else {//not dc, straight up center
                        return [Math.PI / 2, 3 * Math.PI / 2];
                    }
                } else {//down arrow
                    return [4 * Math.PI / 3, 2 * Math.PI / 3];
                }
            } else {//marker to pheno
                return [0, Math.PI];
            }
        }
        
        //create new scale for arrow width in the mini-network
        var netMinEff = Math.min(minEff, minEffBar)
        var netMaxEff = Math.min(maxEff, maxEffBar)
        var phenoArrWidth = d3.scale.linear()
            .domain([netMinEff, netMaxEff])
            .range([1, .005 * svgSize])
        ;
        
        var netDiags = netSvg.selectAll(".netDiag") //draw arrows connecting the nodes
            .data(mAndP)
          .enter().append("path")
            .attr("class", "netDiag")
            .attr("d", function(d, i) {//cubic bezier with an arrowhead
                var angs = calcAngle(mAndP[i].source.name, mAndP[i].target.name, mAndP[i].dc);
                var source = getOnEll(mAndP[i].source.name, angs[0], 0, 0);
                var target = getOnEll(mAndP[i].target.name, angs[1], 0, 0);
                var bez1 = [(source[0] + target[0]) / 2, source[1]];
                var bez2 = [(source[0] + target[0]) / 2, target[1]];
                if(mAndP[i].dc) {//doubly connected, heads don't point out of ellipse at the angle that the line touches
                    if(angs[1] > Math.PI) {//arrow head points down
                        var dire = 1;
                    } else {//arrow head points up
                        var dire = -1;
                    }
                    var tri1 = [target[0] - 5, target[1] + dire * 5];
                    var tri2 = [target[0] + 5, target[1] + dire * 5];
                } else {//head points out of ellipse radially
                    var tri1 = getOnEll(mAndP[i].target.name, angs[1], 7, .07);
                    var tri2 = getOnEll(mAndP[i].target.name, angs[1], 7, -.07);
                }
                var out = 
                    "M" + source[0] + "," + source[1] +
                    "C" +   bez1[0] + "," +   bez1[1] +
                    " " +   bez2[0] + "," +   bez2[1] +
                    " " + target[0] + "," + target[1] +//cubic bezier finishes
                    "L" +   tri1[0] + "," +   tri1[1] +//arrowhead begins
                    "M" + target[0] + "," + target[1] +
                    "L" +   tri2[0] + "," +   tri2[1]
                ;
                return out;
            })
            .attr("fill", "none")
            .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
            .attr("stroke", greenOrRed)
            .attr("stroke-width", function(d, i) {return phenoArrWidth(Math.abs(d.eff));})
            .append("svg:title")
            .text(function(d, i) {
                return "source: " + d.source.name + "\ntarget: " + d.target.name + "\neff: " + d.eff + "\nsig: " + d.sig + "\nqval: " + d.qval;
            })
        ;
    }
    /////////////////////////////END FUNCTIONS//////////////////////////////////////////////
    
    /////////////////BEGIN DEFINITIONS AND DATA FORMATTING//////////////////////////////////
    var nCovar     = capeOut.nCovar;
    var links      = capeOut.data;
    var phenos     = capeOut.phenos;
    var chromoLens = capeOut.chromoLens;
    var pVal       = capeOut.pVal;
    
    d3.select("#pVal")
        .attr("value", pVal)
        .attr("max", pVal)
    ;
    
    var ipv = document.getElementById("pVal")
    ipv.onchange = function(ev) {
        console.log(ipv.value);
    }
    
    //check if touch screen device
    var touch = false;
    if('ontouchstart' in document) {
        touch = true;
    }
    
    var ops = [];     //stores opacities so they can be brought back in sharpOthers
    var hold = false; //start out holding nothing
    var holdOn = -1;  //-1 means holding nothing
    
    var arrowPlots = false;
    var arrow = -1; //none
    
    var svgSize = d3.min([window.innerWidth * .49, window.innerHeight * .99]);
    var textSize = .02 * svgSize;
    var outerRad = .3 * svgSize;
    var innerRad = .29 * svgSize;
    var chromoPad = -.025;        //padding between chromosomes in radians. negative because clockwise
    var triSize = .01 * svgSize;  //arrowhead size
    var triWidth = .015;          //arrowhead width
    
    var vChartPad = .08 * svgSize; //vertical distance beneath each chart
    var hChartPad = .08 * svgSize; //horizontal distance on left of each chart
    
    var phenoDat = [];
    var covars = [];
    var toRemove = [];
   
    //find max effect size and max sig
    var maxEff = 0;
    var maxSig = 0;
    var minSig = Infinity;
    for (var i = 0; i < links.length; i++) {
        //max sig
        if (links[i].sig < minSig) {
            minSig = links[i].sig;
        }
        if (links[i].sig > maxSig) {
            maxSig = links[i].sig;
        }
    }
    
    var totalBP = d3.sum(chromoLens); //calc total base pairs
    var theta = d3.scale.linear()     //linear scale for angles
        .domain([0, totalBP])         //in radians
        .range([0, 2 * Math.PI])      //clock wise
    ;
    
    //find doubly connected markers and note them. This is O(n^2) time
    for (var i = 0; i < links.length; i++) {
        for(var j = 0; j < i; j++) {//check i against all links before it
            if(links[i].source.name == links[j].target.name && links[i].target.name == links[j].source.name) {
                links[j].dc = -1; //dc stands for doubly connected
                links[i].dc =  1; //1 and -1 for add and subtract an offset
            }
        }
    }
    
     //create coordinates for markers and bezier points
    for (var i = 0; i < links.length; i++) {
        //input coordinates for markers
        var sourceCoors = markerPos(links[i].source);
        links[i].source.r = sourceCoors[0];
        links[i].source.a = sourceCoors[1];
        var targetCoors = markerPos(links[i].target);
        links[i].target.r = targetCoors[0];
        links[i].target.a = targetCoors[1];
        
        //create bezier points.
        if(links[i].dc) {//split the arrows of doubly connected markers
            links[i].bez1 = {r: (.05 * links[i].dc + .5) * innerRad, a: links[i].source.a};
            links[i].bez2 = {r: (.05 * links[i].dc + .5) * innerRad, a: links[i].target.a};
        } else {//not doubly connected
            links[i].bez1 = {r: .5 * innerRad, a: links[i].source.a};
            links[i].bez2 = {r: .5 * innerRad, a: links[i].target.a};
        }
    }

    //det max eff btw markers and seperate geno and pheno data
    var minEff = Infinity;
    for(var i = 0; i < links.length; i++) {
        if (Math.abs(links[i].eff) < minEff) {
            minEff = Math.abs(links[i].eff);
        }
        if (Math.abs(links[i].eff) > maxEff) {
            maxEff = Math.abs(links[i].eff);
        }
        if (isNaN(links[i].target.r)) { 
            phenoDat.push(links[i]); //add to pheno
            toRemove.push(i); //will remove from geno
        }
    }
    
    //remove pheno from geno. traverse backwards to avoid indexing issues.
    for(var i = toRemove.length - 1; i >= 0 ; i--) {
        links.splice(toRemove[i], 1);
    }
    
    //count and ID phenotypes and det max eff on pheno
    var maxEffBar = 0;
    var minEffBar = Infinity;
    for(var i = 0; i < phenoDat.length; i++) {
        if (Math.abs(phenoDat[i].eff) > maxEffBar) {
            maxEffBar = Math.abs(phenoDat[i].eff);
        }
        if (Math.abs(phenoDat[i].eff) < minEffBar) {
            minEffBar = Math.abs(phenoDat[i].eff);
        }
    }
    
    var barMax = (.199 * svgSize) / phenos.length; //compute max bar size
    
    var effToWid = d3.scale.linear() //effect to width of arrow
        .domain([minEff, maxEff])
        .range([1, .005 * svgSize])
    ;
        
    var sigToAlp = d3.scale.linear() //significance to alpha (opacity)
        .domain([minSig, maxSig])
        .range([.2, 1])
    ;
    var downPos = [];
    var scale = 1;

    /*var fisheye = d3.fisheye.circular()
        .radius(100)
        .distortion(10)
    ;*/
    ///////////////////END DEFINITIONS AND DATA FORMATTING/////////////////////////////
    
    //////////////////////////////BEGIN DRAWING////////////////////////////////////////
    var svg = d3.select("#iris")
      /*.append("div")
        .attr("id", "overlay")
        .attr("width", svgSize)
        .attr("height", svgSize)
        .attr("pointer-events", "all")
        .call(d3.behavior.zoom().scaleExtent([1, 5]).on("zoom", reZoom))
        .on("dblclick.zoom", null)//don't zoom on double click*/
      .append("svg")   //create svg for circular chromosomes
        .attr("id", "irisSvg")
        .attr("width", svgSize)
        .attr("height", svgSize)
        .call(d3.behavior.zoom().scaleExtent([1, 5]).on("zoom", reZoom))
        .on("dblclick.zoom", null)//don't zoom on double click
        /*.on("mousemove", function() { //fisheye
            fisheye.focus(d3.mouse(document.getElementById("translate")));
            d3.selectAll(".arrows").attr("d", function(d, i) {//cubic bezier curve. 4 points
                bez1 = fisheye({x: d.bez1.r * Math.cos(d.bez1.a), y: d.bez1.r * Math.sin(d.bez1.a)})
                bez2 = fisheye({x: d.bez2.r * Math.cos(d.bez2.a), y: d.bez2.r * Math.sin(d.bez2.a)})
                var out =  
                "M" + d.source.r * Math.cos(d.source.a) + " " + d.source.r * Math.sin(d.source.a) +
                "C" + bez1.x  + " " + bez1.y +
                " " + bez2.x  + " " + bez1.y +
                " " + (d.target.r - triSize) * Math.cos(d.target.a) + " " + (d.target.r - triSize) * Math.sin(d.target.a);
                return out;
            })
            svg.attr("transform", function(d) {
                var fe = fisheye(d); //{x: ..., y: ...
                return "translate(" + [fe.x, fe.y] + ")" + "scale(" + fe.z + ")";
            })
        })*/

      .append("g")
        .attr("id", "zoom")
      .append("g")                    //center it and rotate
        .attr("id", "translate")
        .attr("transform", "translate(" + svgSize * .5 + "," + svgSize * .5 + ")")
      .append("g")
        .attr("id", "rotate")
        .attr("transform", "rotate(" + 90 + ")")
    ;
    d3.select("svg").node().oncontextmenu = function(){return false;}; //turn off right-click menu

    function reZoom() {
        if(d3.event.scale != scale && !touch) { //rescale arrows when zoom level changes
            scale = d3.event.scale;
            //triWidth = triWidth / scale;
            svg.selectAll(".arrows").attr("stroke-width", function(d, i) {
                if(scale < 3) {
                    return effToWid(Math.abs(d.eff)) / scale;
                } else {
                    return effToWid(Math.abs(d.eff)) / 3;
                }
            })
            /*d3.selectAll(".arrowHeads").attr("points", function(d, i) {
                var out = 
                d.target.r * Math.cos(d.target.a) + "," + d.target.r * Math.sin(d.target.a) + " " +
                (d.target.r - triSize) * Math.cos(d.target.a - triWidth) + "," + (d.target.r - triSize) * Math.sin(d.target.a - triWidth) + " " +
                (d.target.r - triSize) * Math.cos(d.target.a + triWidth) + "," + (d.target.r - triSize) * Math.sin(d.target.a + triWidth);
                return out;
            })*/
        }
        d3.select("#zoom").attr("transform", "translate(" + d3.event.translate + ")" + " scale(" + d3.event.scale + ")");
    }
    
    var arrows = svg.selectAll("path.arrows") //arrows
        .data(links)
      .enter().append("path")
        .attr("class", "arrows")
        .attr("id", function(d, i) {return "arrow0" + i;})
        .attr("d", function(d, i) {//cubic bezier curve. 4 points
            var out =  
            "M" + d.source.r * Math.cos(d.source.a) + " " + d.source.r * Math.sin(d.source.a) +
            "C" +   d.bez1.r * Math.cos(d.bez1.a)   + " " +   d.bez1.r * Math.sin(d.bez1.a)   + 
            " " +   d.bez2.r * Math.cos(d.bez2.a)   + " " +   d.bez2.r * Math.sin(d.bez2.a)   +
            " " + (d.target.r - triSize) * Math.cos(d.target.a) + " " + (d.target.r - triSize) * Math.sin(d.target.a);
            return out;
        })
        .attr("stroke", greenOrRed)
        .attr("stroke-width", function(d, i) {return effToWid(Math.abs(d.eff));})
        .attr("fill", "none")
        .attr("transform", "rotate(" + -90 + ")") //not sure why this is necessary, but it makes it work.
        .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
        .on("mousedown", function(d, i) {
            downPos = d3.mouse(document.getElementById("irisSvg"));
        }) //save position of mouse down. used for determining if a drag occured
        .on("mouseup", toggleArrowPlots) //turn on or off the plots connected to that arrow
      .append("svg:title")
        .text(function(d, i) {
            return "source: " + d.source.name + "\ntarget: " + d.target.name + "\neff: " + d.eff + "\nsig: " + d.sig + "\nqval: " + d.qval;
        })
    ;
        
    var arrowHeads = svg.selectAll("polygon.arrowHeads") //arrow heads
        .data(links)
      .enter().append("polygon")
        .attr("class", "arrowHeads")
        .attr("id", function(d, i) {return "arrowHead0" + i})
        .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
        .attr("points", function(d, i) {
            var out = 
            d.target.r * Math.cos(d.target.a) + "," + d.target.r * Math.sin(d.target.a) + " " +
            (d.target.r - triSize) * Math.cos(d.target.a - triWidth) + "," + (d.target.r - triSize) * Math.sin(d.target.a - triWidth) + " " +
            (d.target.r - triSize) * Math.cos(d.target.a + triWidth) + "," + (d.target.r - triSize) * Math.sin(d.target.a + triWidth);
            return out;
        })
        .attr("fill", greenOrRed)
        .attr("transform", "rotate(" + -90 + ")")
    ;
        
    var chromoArc = d3.svg.arc() //create specs for arcs
        .outerRadius(innerRad)
        .innerRadius(outerRad)//backwards so that labels can be put inside
        .startAngle(function(d, i) {return theta(sumTo(chromoLens, i));})
        .endAngle(function(d, i) {return theta(sumTo(chromoLens, i + 1)) + chromoPad;})
    ;    
        
    var chromoArcs = svg.selectAll("path.chromoArcs") //put arcs in svg in form of path
        .data(chromoLens)
      .enter().append("path")     
        .attr("class", "chromoArcs")
        .attr("id", function(d, i) {return "chromoArc0" + i;})
        .attr("d", chromoArc)
        .attr("stroke", "black")
        .attr("fill", "black")
        .attr("opacity", 1)
    ;

    var barHeight = d3.scale.linear() //scale for phenotype bars
        .domain([minEffBar, maxEffBar])
        .range([4, barMax])
    ;
    
    var phenoBar = d3.svg.arc() //radial pheno bars arc specs
        .outerRadius(function(d, i) {return phenos.indexOf(d.target.name) * barMax + outerRad + barHeight(Math.abs(d.eff));})
        .innerRadius(function(d, i) {return phenos.indexOf(d.target.name) * barMax + outerRad;})
        .startAngle(function(d, i) {return d.source.a - .007;})
        .endAngle(function(d, i) {return d.source.a + .007;})
    ;
        
    var phenoBars = svg.selectAll("path.phenoBars") //bars outside of chromosomes showing effect on phenotype as svg path
        .data(phenoDat)
      .enter().append("path")
        .attr("class", "phenoBars")
        .attr("id", function(d, i) {return "phenoBar0" + d.source.chromosome;})
        .attr("d", phenoBar)
        .attr("fill", greenOrRed)
        .attr("stroke", "white")
        .attr("stroke-width", .5)
        .attr("opacity", function(d, i) {return sigToAlp(d.sig);})
    ;
    
    for(var p = 0; p < phenos.length + 1; p++) { //rings to seperate phenoBars by phenotype
        var zeroAxis = d3.svg.arc()
            .outerRadius(function(d, i) {return p * barMax + outerRad;})
            .innerRadius(function(d, i) {return p * barMax + outerRad - .1;})
            .startAngle(function(d, i) {return theta(sumTo(chromoLens, i));})
            .endAngle(function(d, i) {return theta(sumTo(chromoLens, i + 1)) + chromoPad;})
        ;
            
        var zeroAxes = svg.selectAll("path.zeroAxes" + p)
            .data(chromoLens)
          .enter().append("path")
            .attr("class", function(d, i) {return "zeroAxes" + p;})
            .attr("id", function(d, i) {return "zeroAxis" + p +"0" + i;})
            .attr("d", zeroAxis)
            .attr("stroke", "black")
            .attr("opacity", 1)//needs to have an opacity, so it is known what to return to in sharpOthers
            .attr("fill", "black")
        ;
        
        /*if(p % 2 != 0 && p != phenos.length && phenos.length >= 3) { //p is odd. color these rings darker. zebra striping. can't zebra stripe, it inflates opacities
            var zebraStripe = d3.svg.arc()
                .outerRadius((p + 1) * barMax + outerRad)
                .innerRadius(p * barMax + outerRad - .1)
                .startAngle(0)
                .endAngle(2 * Math.PI)
            ;
            
            svg.append("path")
                .attr("class", "zebraStripe")
                .attr("d", zebraStripe)
                .attr("fill", "black")
                .attr("stroke", "none")
                .attr("opacity", .05)
            ;
        }*/
    }
    
    var chromoLabels = svg.selectAll("text.chromoLabels") //label the chromosomes
        .data(chromoLens)
      .enter().append("text")
        .attr("class", "chromoLabels")
        .attr("id", function(d, i) {return "chromoLabel0" + i;})
        .attr("font-size", outerRad - innerRad + "pt")
        .attr("fill", "white")
      .append("textPath")
        .attr("startOffset", "22%")//not 25 because arcs have a width. 22 is a decent guess.
        .attr("text-anchor", "middle")
        .attr("xlink:href", function(d, i) {return "#chromoArc0" + i;})
        .text(function(d, i) {
            if(i >= chromoLens.length - nCovar)//is covar
                return ""; //not enough room
            else {
                return i + 1;
            }
        })
    ;
    
    var phenoLabelArc = d3.svg.arc()
        .innerRadius(function(d, i) {return i * barMax + outerRad + barMax / 3;})
        .outerRadius(function(d, i) {return i * barMax + outerRad + barMax / 3;})
        .startAngle(5 * Math.PI / 4)
        .endAngle(7 * Math.PI / 4)
    ;
    
    var phenoLabelArcs = svg.selectAll("path.phenoLabelArc")
        .data(phenos)
      .enter().append("path")
        .attr("class", "phenoLabelArc")
        .attr("id", function(d, i) {return "phenoLabelArc0" + i;})
        .attr("d", phenoLabelArc)
        .attr("fill", "none")
        .attr("stroke", "none")
    ;
    
    var phenoLabels = svg.selectAll("text.phenoLabels") //label the phenotypes
        .data(phenos)
      .enter().append("text")
        .attr("class", "phenoLabels")
        .attr("opacity", .7)
        .attr("font-size", textSize * 1.75)
        .attr("id", function(d, i) {return "phenoLabel0" + i;})
      .append("textPath") //put label on top. There is nowhere that is guaranteed to be data free.
        .attr("startOffset", "25%")
        .attr("text-anchor", "middle")
        .attr("xlink:href", function(d, i) {return "#phenoLabelArc0" + i;})
        .text(function(d, i) {return d;})
    ;
        
    var mouseRegion = d3.svg.arc() //one for each chromosome. for clicking and mouseover interaction.
        .outerRadius(outerRad + phenos.length * barMax)
        .innerRadius(innerRad)
        .startAngle(function(d, i) {return theta(sumTo(chromoLens, i)) + chromoPad / 2;})
        .endAngle(function(d, i) {return theta(sumTo(chromoLens, i + 1)) + chromoPad / 2;})
    ;
    
    var chromoDown = [];
    var mouseRegions = svg.selectAll("path.mouseRegions") //chromosomal regions for mouse interaction
        .data(chromoLens)
      .enter().append("path")
        .attr("id", function(d, i) {return "mouseRegion0" + i;})
        .attr("class", "mouseRegions")
        .attr("d", mouseRegion)
        .attr("fill", "white")
        .attr("opacity", .2)
        .attr("stroke", "none")
        .on("mouseover", fadeOthers) //fade objects not connected to the moused over one
        .on("mouseout", sharpOthers) //bring them back
        .on("mousedown", function(d, i) {
            chromoDown = d3.mouse(document.getElementById("irisSvg"));
        })
        .on("mouseup", toggleHold)     //toggle hold of current opacity state
    ;
    if(touch) {
        d3.selectAll("path.mouseRegions")
            .on("mouseover", null)
            .on("mouseout", null)
        ;
    }
    ///////////////////////////////END DRAWING/////////////////////////////////////////////////////
}
