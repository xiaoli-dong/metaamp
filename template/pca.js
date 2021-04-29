


var charts = [];
var getChartConfig = function(renderId, titletxt, data, type, xtitle, ytitle, legend) {
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

        plotBackgroundColor: 'rgba(255, 255, 255, .9)'

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
        plotLines: [{
            color: '#0000ff',
            width: 2,
            value: 0,
            id: 'plotline-1'
        }],
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
        plotLines: [{
            color: '#0000ff',
            width: 2,
            value: 0,
            id: 'plotline-1'
        }],
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


    config.series = data;



    return config;
};

var formatNumber = d3.format(",.2f"); // zero decimal places
format = function(d) {
    return formatNumber(d);
};


$(document).ready(function() {
    axes = axesdata;
    loading = loadingdata;
    var p12 = axes.map(function(d) {
        return {
            sname: d.group,
            x: +d.axis1,
            y: +d.axis2,
        };
    });

    console.log(p12);

    var p13 = axes.map(function(d) {
        return {
            sname: d.group,
            x: +d.axis1,
            y: +d.axis3,
        };
    });

    console.log(p13);

    var p23 = axes.map(function(d) {
        return {
            sname: d.group,
            x: +d.axis2,
            y: +d.axis3,
        };
    });

    var p1_loading = loading[0].loading;
    var p2_loading = loading[1].loading;
    var p3_loading = loading[2].loading;

    var p12data = [{
        name: '',
        color: 'rgba(119, 152, 191, .5)',
        data: p12

    }];

    var p13data = [{
        name: '',
        color: 'rgba(119, 152, 191, .5)',
        data: p13
    }];
    var p23data = [{
        name: '',
        color: 'rgba(119, 152, 191, .5)',
        data: p23
    }];

    var p1_title = 'P1 - percent variation explained ' + format(p1_loading) + '%';
    var p2_title = 'P2 - percent variation explained ' + format(p2_loading) + '%';
    var p3_title = 'P3 - percent variation explained ' + format(p3_loading) + '%';
    //now, creating a new chart is easy!

    var chart1 = new Highcharts.Chart(getChartConfig("pca12", "PCoA P1 vs P2", p12data, "scatter", p1_title, p2_title, false));
    charts.push(chart1);
    var chart2 = new Highcharts.Chart(getChartConfig("pca13", "PCoA P1 vs P3", p13data, "scatter", p1_title, p3_title, false));
    charts.push(chart2);
    var chart3 = new Highcharts.Chart(getChartConfig("pca23", "PCoA P2 vs P3", p23data, "scatter", p2_title, p3_title, false));
    charts.push(chart3);


    // the button action
    $('#button1').click(function() {
        var opt = chart1.series[0].options;
        opt.dataLabels.enabled = !opt.dataLabels.enabled;
        chart1.series[0].update(opt);
    });

    // the button action
    $('#button2').click(function() {
        var opt = chart2.series[0].options;
        opt.dataLabels.enabled = !opt.dataLabels.enabled;
        chart2.series[0].update(opt);
    });

    // the button action
    $('#button3').click(function() {
        var opt = chart3.series[0].options;
        opt.dataLabels.enabled = !opt.dataLabels.enabled;
        chart3.series[0].update(opt);
    });


    



});
