
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


var margin = { top: 100, right: 100, bottom: 100, left: 100 },
cellsize=50,
width = cellsize * ccount,
height = cellsize * rcount;


domain = d3.set(data.map(function(d) {
    return d.colid
})).values(),
num = Math.sqrt(data.length),

color = d3.scaleLinear()
    .domain([-1, 0, 1])
    .range(["#B22222", "white", "#000080"])

var x = d3.scalePoint()
    .range([0, width])
    .domain(domain),

y = d3.scalePoint()
    .range([0, height])
    
    .domain(domain),
xSpace = x.range()[1] - x.range()[0],
ySpace = y.range()[1] - y.range()[0];


var svg = d3.select("#chart")
    .append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")")
    
    
var cor = svg.selectAll(".cor")
    .data(data)
    .enter()
    .append("g")
    .attr("class", "cor")
    .attr("transform", function(d) {
        return "translate(" + x(d.colid) + "," + y(d.rowid) + ")";
    })
    


cor.append("rect")
    .attr("width", xSpace/ccount)
    .attr("height", ySpace/rcount)
    .attr("x", -xSpace / num/2)
    .attr("y", -ySpace / num/2)
    
cor.filter(function(d){
    var ypos = domain.indexOf(d.rowid);
    var xpos = domain.indexOf(d.colid);
    for (var i = (ypos + 1); i < num; i++){
        if (i === xpos) return false;
    }
    return true;
})
    .append("text")
    .attr("y", 5)
    .text(function(d) {
        if (d.colid === d.rowid) {
            return d.colid;
        } else {
            return parseFloat(d.value).toFixed(2);
        }
    })
    .style("fill", function(d){
        
	if(d.value == 0) {
            return "#000";
        }
	else {
            return color(d.value);
        }
    });

cor.filter(function(d){
    
    var ypos = domain.indexOf(d.rowid);
    var xpos = domain.indexOf(d.colid);
    for (var i = (ypos + 1); i < num; i++){
        //if (i == xpos && d.value != 0) return true;
	if (i == xpos) return true;
    }
    return false;
})
    .append("circle")
    .attr("r", function(d){
        //return (width / (num * 2 + 2)) * (Math.abs(d.value) + 0.1);
	return (cellsize/3) * (Math.abs(d.value));
    })
    .style("fill", function(d){
        if (d.value == 0) {
            return "white";
        } else {
            return color(d.value);
        }
    })
.attr("stroke", function(d){
        if (d.value == 0) {
            return "white";
        } else {
	    
	    return "grey";
        }
    })

cor.on("mouseover", function(d){
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


var legend_cell_count = 10;

var margin_legend = { top: 10, right: 10, bottom: 50, left: 100 },
//width_legend = 50 * legend_cell_count,
//width_legend = width,
width_legend = 650,
height_legend = 100;


var svg_legend = d3.select("#svg").append("svg")
    .attr("width", width + margin.left + margin.right)
    .attr("height", height_legend)

var shape_width = Math.floor(width_legend/legend_cell_count);

svg_legend.append("g")
    .attr("class", "legendLinear")
    .attr("transform", "translate(" + margin_legend.left + "," + margin_legend.top + ")");

var legendLinear = d3.legendColor()
    .shapeWidth((shape_width))
    .cells(legend_cell_count )
    .orient('horizontal')
    .scale(color)
    .labelFormat(d3.format(".2f"))

svg_legend.select(".legendLinear")
    .call(legendLinear);



