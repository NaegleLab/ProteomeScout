$(document).ready( function() {
	d3.selectAll(".experiment_data")
		.each(function() { processRuns(this); });
});

function processDataPoint(dp, points, stddev) {
	csv = d3.select(dp).text();
	columns = csv.split(",");
	
	x = columns[0];
	y = columns[1];
	type = columns[2];
	
	if(y=="None"){
		return;
	}
	
	if(type == "stddev"){
		ypos = -1;
		
		for(j = 0; j < points.length; j++) {
			if(points[j].label == x)
				ypos = points[j].y
		}
		if(ypos == 0 && console != undefined)
			console.log("Warning: data point for sigma at {0} not found".format(x));
		
		stddev.push({'label':x, 'y':ypos, 'dev':y})
	} else {
		points.push({'label':x,'y':y})
	}
}

function processRun(run, experiment_data, run_map) {
	name = d3.select(run).select(".name").text();
	
	points = []
	stddev = []
	
	phosphopeps = 
		d3.select(run)
			.select(".phosphopeps")
			.text();
	
	units = 
		d3.select(run)
			.select(".units")
			.text();
	
	d3.select(run)
		.selectAll(".datapoint")
		.each(function() { processDataPoint(this, points, stddev); });
	
	isTime = (units.indexOf("time") >= 0)
	
	if (! (name in run_map)){
		run_map[name] = {'name':name, 'series':[], 'isTime': isTime, 'axis': units};
	}
	
	run_map[name].series.push({'peps':phosphopeps, 'points':points, 'stddev':stddev});
}

function processRuns(experiment_data){
	run_map = {}
	
	d3.select(experiment_data)
		.selectAll(".run")
		.each(function() { processRun(this, experiment_data, run_map); });
	
	createGraphs(experiment_data, run_map);
}

function createGraphs(experiment_data, run_map) {
	for(r in run_map) {
		run = run_map[r]
		if(run.isTime)
			createTimeSeriesGraph(experiment_data, run);
		else
			createBarGraph(experiment_data, run);
	}
}

function getMaxValue(run, property) {
	mval = -10000.0;
	for(i = 0; i < run.series.length; i++) {
		for(j = 0; j < run.series[i].points.length; j++) {
			val = parseFloat(run.series[i].points[j][property]);
			if(val > mval) {
				mval = val;
			}
		}
	}
	return mval;
}

function createTimeSeriesGraph(experiment_data, run) {
	var w = 475;
	var h = 400;
	var margin = 40;
	var rmargin = 130;
	var ceiling = 1.15;
	var lwidth = 125;
	var colors = d3.scale.category20();
	
	container = d3.select(experiment_data).select(".chart").append("span");
	parent = container.append("span");
	graph = createGraph(parent, run.name, w, h, margin);
	
	xmax = getMaxValue(run, "label");
	ymax = getMaxValue(run, "y") * ceiling;
	
	xaxis = d3.scale.linear().domain([0, xmax]).range([margin, w-rmargin]);
	yaxis = d3.scale.linear().domain([0, ymax]).range([h-margin, margin]);
	
	legendEntries = [];
	
	xticks = []
	for(i = 0; i < run.series.length; i++){
		pts = []
		stddev = []
		for(j = 0; j < run.series[i].points.length; j++){
			pt = run.series[i].points[j];
			
			xticks.push( parseFloat(pt.label) );
			pts.push( { 'x':parseFloat(pt.label), 'y':parseFloat(pt.y) } );
		}
		
		for(j = 0; j < run.series[i].stddev.length; j++){
			pt = run.series[i].stddev[j]
			
			stddev.push( {'x': parseFloat(pt.label), 'y': parseFloat(pt.y), 's': parseFloat(pt.dev)} )
		}
		name = run.series[i].peps
		
		addTimeSeries(graph, name, pts, xaxis, yaxis, colors(i));
		addErrorBars(graph, name, stddev, xaxis, yaxis, colors(i));
		
		legendEntries.push({'name':name, 'color':colors(i)});
	}
	
	addAxes(graph, run.axis, Array.unique(xticks), yaxis.ticks(7), xaxis, yaxis, false);
	addLegend(graph, legendEntries, w-lwidth, margin, lwidth, "line");
	addExport(parent, container);
}

function getArray(run, property) {
	rval = [];
	for(i = 0; i < run.series.length; i++) {
		for(j = 0; j < run.series[i].points.length; j++) {
			rval.push(run.series[i].points[j][property]);
		}
	}
	return rval;
}

function createBarGraph(experiment_data, run) {
	var h = 440;
	var margin = 40;
	var rmargin = 130;
	var bmargin = 80;
	var ceiling = 1.15;
	var lwidth = 125;
	var colors = d3.scale.category20();

	xvals = Array.unique(getArray(run, "label"));
	xvals.unshift("");
	xvals.push(" ");
	
	defaultBarWidth = 30;
	
	var w = margin + rmargin + defaultBarWidth * run.series.length * xvals.length;
	
	if(w > 750)
		w = 750;
	
	parent = d3.select(experiment_data).select(".chart");
	graph = createGraph(parent, run.name, w, h, margin);
	
	
	ymax = getMaxValue(run, "y") * ceiling;
	
	xaxis = d3.scale.ordinal().domain(xvals).rangePoints([margin, w-rmargin]);
	yaxis = d3.scale.linear().domain([0, ymax]).range([h-bmargin, margin]);
	
	legendEntries = [];
	
	xticks = []
	for(i = 0; i < run.series.length; i++){
		pts = [];
		stddev = [];
		
		for(j = 0; j < run.series[i].points.length; j++){
			pt = run.series[i].points[j];
			xticks.push( parseFloat(pt.label) );
			pts.push( { 'x':pt.label, 'y':parseFloat(pt.y) } );
		}
		
		for(j = 0; j < run.series[i].stddev.length; j++){
			pt = run.series[i].stddev[j]
			stddev.push( {'x':pt.label, 'y': parseFloat(pt.y), 's': parseFloat(pt.dev)} )
		}
		name = run.series[i].peps
		
		addBarSeries(graph, name, pts, xaxis, yaxis, colors(i), i, run.series.length);
		addErrorBars(graph, name, stddev, xaxis, yaxis, "#000", i, run.series.length);
		
		legendEntries.push({'name':name, 'color':colors(i)});
	}
	
	addAxes(graph, run.axis, xvals, yaxis.ticks(7), xaxis, yaxis, true);
	addLegend(graph, legendEntries, w-lwidth, margin, lwidth, "square");
	addExport(parent, d3.select(experiment_data));
}