
var charts = [];
var getChartConfig = function(renderId, titletxt, mydata, type, xtitle, ytitle, legend, stress, rsquare) {
    var config = {};
    config.colors = ['#f45b5b', '#8085e9', '#8d4654', '#7798BF', '#aaeeee', '#ff0066', '#eeaaee',
		     '#55BF3B', '#DF5353', '#7798BF', '#aaeeee'],

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
    config.subtitle = {
	text: 'stress: ' + stress + ', R-squared: ' + rsquare,
	margin: 30,

        style: {
            color: '#666',
            font: 'bold 12px "Trebuchet MS", Verdana, sans-serif'
        }
    };
    config.xAxis = {
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
            enabled: true,
            text: xtitle,
            style: {
                color: '#333',
                fontWeight: 'bold',
                fontSize: '12px',
                fontFamily: 'Trebuchet MS, Verdana, sans-serif'

            }
        }

    };

    config.yAxis = {
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
    };
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
                formatter: function() {
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
        formatter: function() {
            return '<strong>' + this.point.sname + '</strong>';
        },
        style: {
            color: '#666666',
            fontWeight: 'bold',
            font: '11px "Trebuchet MS", Verdana, sans-serif'
        }

    };


    config.series = [
	{
	    name: '',
	    data:mydata
	}
    ];



    return config;
};

var formatNumber = d3.format(",.2f"); // zero decimal places
format = function(d) {
    return formatNumber(d);
};

$(document).ready(function() {
    

    stress = stressdata;
    axes = axesdata;
    
    
    var str_and_rsq = stress.map(function(d) {
        return {
            Stress: +d.Stress,
            Rsq: +d.Rsq
        };
    });
    var sorted_stress =  str_and_rsq.sort(function(x, y){
	
	return d3.ascending(x.index, y.index);
    });
    var stress_value = sorted_stress[0].Stress;
    var rsq = sorted_stress[0].Rsq;

    var axis12 = axes.map(function(d) {
        return {
            sname: d.group,
            x: +d.axis1,
            y: +d.axis2
        };
    });



    console.log(axis12);

    var chart1 = new Highcharts.Chart(getChartConfig("nmds", "NMDS ordination", axis12, "scatter", "Axis1", "Axis2", false, format(stress_value), format(rsq)));
    charts.push(chart1);
    // the button action
    $('#button1').click(function() {
        var opt = chart1.series[0].options;
        opt.dataLabels.enabled = !opt.dataLabels.enabled;
        chart1.series[0].update(opt);
    });
    

});


