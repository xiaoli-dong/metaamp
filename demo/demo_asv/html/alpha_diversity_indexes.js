

var charts = [];
var getChartConfig = function (renderId, sname, titletxt, mydata, dist,
			       cat, type, xtitle, ytitle, legend, linked, linked_data) {
    var config = {};
    config.chart = {
	renderTo: renderId,
	defaultSeriesType: type,
	borderWidth: 2,
	spacingLeft: 20,
	spacingRight: 20,
	spacingTop: 20,
	spacingBottom: 20,
	backgroundColor: {
	    linearGradient: {
		x1: 0,
		y1: 0,
		x2: 1,
		y2: 1
	    },
	    stops: [
		[0, 'rgb(255, 255, 255)'],
		[1, 'rgb(240, 240, 255)']
	    ]
	},
	plotBackgroundColor: 'rgba(255, 255, 255, .9)',
	plotShadow: true,
	plotBorderWidth: 1
    };
    config.credits = {
	enabled: false
    };
    config.title = {
	text: titletxt,
	margin: 30,
	style: {
	    color: '#000',
	    font: 'bold 16px "Trebuchet MS", Verdana, sans-serif'
	}
    };
    config.xAxis = {
	categories: cat,
	gridLineWidth: 1,
	lineColor: '#000',
	tickColor: '#000',
	labels: {
	    style: {
		color: '#000',
		font: '11px Trebuchet MS, Verdana, sans-serif'
	    }
	},
	title: {
	    text: xtitle,
	    style: {
		color: '#333',
		fontWeight: 'bold',
		fontSize: '12px',
		fontFamily: 'Trebuchet MS, Verdana, sans-serif'
	    }
	}
    };
    config.yAxis = [{
	minorTickInterval: 'auto',
	lineColor: '#000',
	lineWidth: 1,
	tickWidth: 1,
	tickColor: '#000',
	labels: {
	    style: {
		color: '#000',
		font: '11px Trebuchet MS, Verdana, sans-serif'
	    }
	},
	title: {
	    text: ytitle,
	    style: {
		color: '#333',
		fontWeight: 'bold',
		fontSize: '12px',
		fontFamily: 'Trebuchet MS, Verdana, sans-serif'
	    }
	}
    }];
    config.legend = {
	enabled: legend
    };
    config.plotOptions = {
	scatter: {
	    marker: {
		radius: 3,
		fillColor: '#ff0000',
		states: {
		    hover: {
			enabled: true,
			lineColor: 'rgb(100,100,100)'
		    }
		}
	    },
	    dataLabels: {
		enabled: true,
		allowOverlap: true,
		formatter: function () {
		    return this.point.sname;
		},
		style: {
		    color: '#666666',
		    fontWeight: 'bold',
		    font: '11px "Trebuchet MS", Verdana, sans-serif'
		}
	    },
	    
	   
	    states: {
		hover: {
		    marker: {
			enabled: false
		    }
		}
	    }
	}
    };
    config.navigation = {
	buttonOptions: {
	    align: 'right'
	}
    };
    
    config.tooltip = {
	style: {
	    color: '#666666',
	    fontWeight: 'bold',
	    font: '11px "Trebuchet MS", Verdana, sans-serif'
	}
    };
	

    if(linked){
	config.series = [
	    {
		name: sname,
		type: 'scatter',		 
		labels: {
		    format: '{value}'
		},
		data: linked_data
		},
		{
		name: sname,
	    labels: {
		format: '{value}'
	    },
	    data: mydata
	},
		     
			];
    }
    else{
	
	config.series = [{
	    name: sname,
	    labels: {
		format: '{value}'
	    },
	    data: mydata
	}	     
			 
			 
			];
    }
    return config;
};
var formatNumber = d3.format(",.2f"); // zero decimal places
format = function (d) {
    return formatNumber(d);
};

$(document)
    .ready(function () {
	
	var cat = data.map(function (d) {
	    return d.group;
	});
	var sob_data = data.map(function (d) {
	    return +d.sobs
	});
	var chao_data_interval = data.map(function (d) {
	    return [+d.chao_lci, +d.chao_hci];
	});
	var chao_data = data.map(function (d) {
	    return [+d.chao];
	});
	
	console.log(chao_data);
	
	var ace_data_interval = data.map(function (d) {
	    return [+d.ace_lci, +d.ace_hci];
	});
	
	var ace_data = data.map(function (d) {
	    return +d.ace;
	});
	
	var jackknife_data_interval = data.map(function (d) {
	    return [+d.jackknife_lci, +d.jackknife_hci];
	});
	var jackknife_data = data.map(function (d) {
	    return +d.jackknife;
	});
	
	var shannon_data_interval = data.map(function (d) {
	    return [+d.shannon_lci, +d.shannon_hci];
	});
	var shannon_data = data.map(function (d) {
	    return +d.shannon;
	});
	var npshannon_data = data.map(function (d) {
	    return +d.npshannon;
	});
	var simpson_data_interval = data.map(function (d) {
	    return [+d.simpson_lci, +d.simpson_hci];
	});
	var simpson_data = data.map(function (d) {
	    return +d.simpson;
	});
	var dis = 10;
	//now, creating a new chart is easy!
	charts.push(new Highcharts.Chart(getChartConfig(
	    "sob", "sob",
	    "sob distribution", sob_data,
	    dis, cat, "lollipop", "Samples",
	    "number of observed OTUs",
	    false)));

	
	charts.push(new Highcharts.Chart(getChartConfig(
	    "chao", "chao",
	    "chao distribution", chao_data_interval,
	    dis, cat, "errorbar", "Samples",
	    "Chao1 richness estimate",
	    false, true, chao_data)));
	
	
	
	charts.push(new Highcharts.Chart(getChartConfig(
	    "ace", "ace",
	    "ace distribution", ace_data_interval,
	    dis, cat, "errorbar", "Samples",
	    "ACE richness estimate", false, true, ace_data)));
	
	charts.push(new Highcharts.Chart(getChartConfig(
	    "jackknife", "jackknife",
	    "jackknife distribution",
	    jackknife_data_interval, dis, cat,
	    "errorbar", "Samples",
	    "Jackknife richness estimate",
	    false, true, jackknife_data)));
	charts.push(new Highcharts.Chart(getChartConfig(
	    "shannon", "shannon",
	    "shannon distribution",
	    shannon_data_interval, dis, cat, "errorbar",
	    "Samples",
	    "Shannon diversity index ",
	    false, true, shannon_data)));
	charts.push(new Highcharts.Chart(getChartConfig(
	    "npshannon", "npshannon",
	    "npshannon distribution",
	    npshannon_data, dis, cat,
	    "lollipop", "Samples",
	    "Non-parametric Shannon diversity index",
	    false)));
	charts.push(new Highcharts.Chart(getChartConfig(
	    "simpson", "simpson",
	    "Simpson distribution",
	    simpson_data_interval, dis, cat, "errorbar",
	    "Samples",
	    "Simpson diversity index",
	    false, true, simpson_data)));
	
    });
