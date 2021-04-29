Array.prototype.contains = function(v) {
  for (var i = 0; i < this.length; i++) {
    if (this[i] === v) return true;
  }
  return false;
};

Array.prototype.unique = function() {
  var arr = [];
  for (var i = 0; i < this.length; i++) {
    if (!arr.contains(this[i])) {
      arr.push(this[i]);
    }
  }
  return arr;
}

//var duplicates = [1, 3, 4, 2, 1, 2, 3, 8];
//var uniques = duplicates.unique(); // result = [1,3,4,2,8]

//console.log(uniques);
/**
 * The "median" is the "middle" value in the list of numbers.
 *
 * @param {Array} numbers An array of numbers.
 * @return {Number} The calculated median value from the specified numbers.
 */
function median(numbers) {
    // median of [3, 5, 4, 4, 1, 1, 2, 3] = 3
    var median = 0, numsLen = numbers.length;
    numbers.sort();

    if (
        numsLen % 2 === 0 // is even
    ) {
        // average of two middle numbers
        median = (numbers[numsLen / 2 - 1] + numbers[numsLen / 2]) / 2;

    } else { // is odd
        // middle number only
        median = numbers[(numsLen - 1) / 2];
    }

    return median;
}

/**
 * The "mean" is the "average" you're used to, where you add up all the numbers
 * and then divide by the number of numbers.
 *
 * For example, the "mean" of [3, 5, 4, 4, 1, 1, 2, 3] is 2.875.
 *
 * @param {Array} numbers An array of numbers.
 * @return {Number} The calculated average (or mean) value from the specified
 *     numbers.
 */
function mean(numbers) {
    var total = 0, i;
    for (i = 0; i < numbers.length; i += 1) {
        total += numbers[i];
    }
    return total / numbers.length;
}
var rlabels = data.map(a=>a.rowid).unique();
var clabels = data.map(a=>a.colid).unique();
var values = [];
(
data.map(a=>a.value)).forEach(element => {
    if(element != 0){
	values.push(element);
    }
});

values.sort(function(a, b){return a-b});

var data_min = values[0];
var data_max = values[values.length-1];
var data_median = median(values.map(i=>Number(i)));
console.log("min=" + data_min + ", median=" + data_median + ",max=" + data_max);

console.log(values);

ccount = clabels.length;
rcount = rlabels.length;
var n = data.length;

console.log(rlabels);
console.log(data);

//************************** define svg parameter***********************
var margin = { top: 100, right: 100, bottom: 100, left: 100 },
    cellsize=35,
    width = cellsize * ccount,
    height = cellsize * rcount;

// Build X scales and axis:
var x = d3.scaleBand()
    .range([ 0, width ])
    .domain(clabels)
    .padding(0.01);

// Build y scales:
var y = d3.scaleBand()
    .range([0, height])
    .domain(rlabels)
    .padding(0.01);

// Build color scale


var chromaScale = chroma.scale(["tomato", "white", "steelblue"])
    .domain([data_min, data_min, data_max]);
    //.domain([0, 0, 33]);


var linear = d3.scaleLinear()
    .domain([data_min, data_median, data_max])
    .range(["#B22222", "white", "#000080"])

var legend_cell_count = 10;

var margin_legend = { top: 10, right: 10, bottom: 50, left: 100 },
    //width_legend = 50 * legend_cell_count,
width_legend = 650,
height_legend = 100;


var svg_legend = d3.select("#svg").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height_legend)

var shape_width = Math.floor(width_legend/legend_cell_count);
svg_legend.append("g")
    .attr("class", "legendLinear")
    //.attr("transform", "translate(" + margin_svg.left + ", 20)");
    .attr("transform", "translate(" + margin_legend.left + "," + margin_legend.top + ")");
var legendLinear = d3.legendColor()
//.shapeWidth(30)
    .shapeWidth(shape_width)
    .cells(legend_cell_count )
    .orient('horizontal')
    .scale(linear);

svg_legend.select(".legendLinear")
    .call(legendLinear);






var rowSortOrder=false;
var colSortOrder=false;
var rowid_SortOrder = false;
var colid_SortOrder = false;

var svg = d3.select("#chart").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");
    //.attr("transform", "translate(20,100)");



var rowLabels = svg.append("g")
    .selectAll(".rowLabelg")
    .data(rlabels)
    .enter()
    .append("text")
    .text(function (d) { return d; })
    .attr("x", 0)
    .attr("y", function (d, i) { return i * cellsize; })
    .style("text-anchor", "end")
    .attr("transform", "translate(-6," + cellsize / 1.5 + ")")
    .attr("class", function (d,i) { return "rowLabel mono r"+i;} )
    .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
    .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
    .on("click", function(d,i) {rowSortOrder=!rowSortOrder; sortbylabel("r",d,rowSortOrder);});


var colLabels = svg.append("g")
    .selectAll(".colLabelg")
    .data(clabels)
    .enter()
    .append("text")
    .text(function (d) { return d; })
    .attr("x", 0)
    .attr("y", function (d, i) { return i * cellsize; })
    .style("text-anchor", "start")
    .attr("transform", "translate("+cellsize/2 + ",-6) rotate (-90)")
    .attr("class",  function (d,i) { return "colLabel mono c"+i;} )
    .on("mouseover", function(d) {d3.select(this).classed("text-hover",true);})
    .on("mouseout" , function(d) {d3.select(this).classed("text-hover",false);})
    .on("click", function(d,i) {colSortOrder=!colSortOrder;  sortbylabel("c",d,colSortOrder);});

var heatMap = svg.append("g").attr("class","g3")
    .selectAll(".cellg")
    .data(data,function(d){return d.rowid + ":" + d.colid;})
    .enter()
    .append("rect")
    .attr("x", function(d) {return x(d.colid) })
    .attr("y", function(d) { return y(d.rowid) })
    .attr("class", function(d){return "cell cell-border cr_"+ d.rowid +" cc_"+ d.colid })
    .attr("width", cellsize)
    .attr("height", cellsize)
    //.style("fill", function(d) { return linear(d.value); })
    .style("fill", function(d) { if(d.value == 0 ){return "grey";}else{return linear(d.value); }})
    //.attr("fill-opacity",function(d){if(d.value == 0 ){return 0.5;}else{return 1}})
    .on("mouseover", function(d){
	//highlight text
	d3.select(this).classed("cell-hover",true);
	d3.selectAll(".rowLabel").classed("text-highlight",function(r,ri){ return ri== d.rowid;});
	d3.selectAll(".colLabel").classed("text-highlight",function(c,ci){ return ci== d.colid});

	//Update the tooltip position and value
	d3.select("#tooltip")
	    .style("left", (d3.event.pageX+10) + "px")
	    .style("top", (d3.event.pageY-10) + "px")
	    .select("#value")
	    .text(d.rowid+","+d.colid + ": " + d.value);
	//Show the tooltip
	d3.select("#tooltip").classed("hidden", false);
    })
    .on("mouseout", function(){
	d3.select(this).classed("cell-hover",false);
	d3.selectAll(".rowLabel").classed("text-highlight",false);
	d3.selectAll(".colLabel").classed("text-highlight",false);
	d3.select("#tooltip").classed("hidden", true);
    })
;



// Change ordering of cells
function sortbylabel(rORc,element,sortOrder){

    var t = svg.transition().duration(1000);
    var values=[];
    var sorted; // sorted is zero-based index
    //alert(".c"+rORc+"_"+element);
    d3.selectAll(".c"+rORc+"_"+element)
	.filter(function(ce){
	    console.log(ce);
	    values.push(ce.value);
	});

    //sort by value in a row
    if(rORc=="r"){
	sorted=d3.range(ccount).sort(function(a,b){ if(sortOrder){ return values[b]-values[a];}else{ return values[a]-values[b];}});
	console.log(sorted);
	var clabels_sorted = [];
	console.log(element);
	sorted.forEach(function(i) {
	    clabels_sorted.push(clabels[i]);
	});

	// Build X scales:
	var xsorted = d3.scaleBand()
	    .range([ 0, width ])
	    .domain(clabels_sorted)
	    .padding(0.01);
	t.selectAll(".cell")
	    .attr("x", function(d, i) { return xsorted(d.colid) });
	t.selectAll(".colLabel")
	    .attr("y", function (d, i) { return sorted.indexOf(i) * cellsize; });
    }

    //sort by coumn
    else{
	sorted=d3.range(rcount).sort(function(a,b){if(sortOrder){ return values[b]-values[a];}else{ return values[a]-values[b];}});
	console.log(sorted);
	var rlabels_sorted = [];
	console.log(element);
	sorted.forEach(function(i) {
	    rlabels_sorted.push(rlabels[i]);
	});

	// Build y scales and axis:
	var ysorted = d3.scaleBand()
	    .range([0, height])
	    .domain(rlabels_sorted)
	    .padding(0.01);

	t.selectAll(".cell")
	    .attr("y", function(d) { return ysorted(d.rowid)});
	t.selectAll(".rowLabel")
	    .attr("y", function (d, i) { return sorted.indexOf(i) * cellsize; });
    }
}


var rowid_SortOrder = true;
var colid_SortOrder = true;

d3.selectAll("input[name='sort']").on("change",function(){

    if(this.value === "rowid"){
	order(this.value, rowid_SortOrder, colid_SortOrder);
	rowid_SortOrder = !rowid_SortOrder;
    }
    else if(this.value === "colid"){
	order(this.value, rowid_SortOrder, colid_SortOrder);
	colid_SortOrder = !colid_SortOrder;
    }
});


function order(value, rSortOrder, cSortOrder){

    var t = svg.transition().duration(1000);

    if (value=="rowid"){

	rlabels_sorted_index = d3.range(rcount);
	if(rSortOrder){
	    rlabels_sorted_index.reverse();
	}
	else{
	    rlabels_sorted_index.sort();
	}


	var rlabels_sorted = [];
	rlabels_sorted_index.forEach(function(i) {
	    rlabels_sorted.push(rlabels[i]);
	});
	console.log("rlabels_sorted=");
	console.log(rlabels_sorted);
	// Build y scales and axis:
	var ysorted = d3.scaleBand()
	    .range([0, height])
	    .domain(rlabels_sorted)
	    .padding(0.01);

	t.selectAll(".cell")
	    .attr("y", function(d) { return ysorted(d.rowid)});
	t.selectAll(".rowLabel")
	    .attr("y", function (d, i) { return rlabels_sorted_index.indexOf(i) * cellsize; });

    }

    else if(value=="colid"){

	clabels_sorted_index = d3.range(ccount);
	if(cSortOrder){
	    clabels_sorted_index.reverse();
	}
	else{
	    clabels_sorted_index.sort();
	}

	console.log(clabels_sorted_index);
	var clabels_sorted = [];
	clabels_sorted_index.forEach(function(i) {
	    clabels_sorted.push(clabels[i]);
	});

	// Build X scales and axis:
	var xsorted = d3.scaleBand()
	    .range([ 0, width ])
	    .domain(clabels_sorted)
	    .padding(0.01);

	t.selectAll(".cell")
	    .attr("x", function(d, i) { return xsorted(d.colid) });
	t.selectAll(".colLabel")
	    .attr("y", function (d, i) { return clabels_sorted_index.indexOf(i) * cellsize; });

    }


}
