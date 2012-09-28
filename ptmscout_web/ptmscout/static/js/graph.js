Array.max = function( array ){
    return Math.max.apply( Math, array );
};
Array.min = function( array ){
    return Math.min.apply( Math, array );
};

Array.unique = function( array ){
	return array.filter(function(itm,i,a){
			return i==a.indexOf(itm);
		});
};

String.prototype.format = function() {
	  var args = arguments;
	  return this.replace(/{(\d+)}/g, function(match, number) { 
	    return typeof args[number] != 'undefined'
	      ? args[number]
	      : match
	    ;
	  });
	};

function addLegend(graph, legendEntries, xpos, ypos, width, marker) {
	h = 20;
	
	height = legendEntries.length * h;
	
	legend = graph.append("g")
				.attr("class", "legend")
				.attr("transform", "translate({0},{1})".format( xpos, ypos ));
	
	legend.append("rect")
		.attr("class", "bg")
		.attr("x", 0)
		.attr("y", 0)
		.attr("width", width)
		.attr("height", height);
	
	legend.selectAll("text.entry")
			.data(legendEntries)
		.enter().append("text")
			.attr("class", "entry")
			.attr("x", 30)
			.attr("y", function(d,i) { return (i+1) * h })
			.attr("dy", -6)
			.text(function(d) { return d.name });
	
	if(marker == "line"){
		legend.selectAll("line.marker")
				.data(legendEntries)
			.enter().append("line")
				.attr("class", "marker")
				.attr("x1", 5)
				.attr("x2", 25)
				.attr("y1", function(d,i) {return (i+1) * h - h/2; })
				.attr("y2", function(d,i) {return (i+1) * h - h/2; })
				.attr("stroke", function(d){ return d.color });
	}
	if(marker == "square"){
		legend.selectAll("rect.marker")
			.data(legendEntries)
		.enter().append("rect")
			.attr("class", "marker")
			.attr("x", 10)
			.attr("width", 15)
			.attr("y", function(d,i) {return i * h+5; })
			.attr("height", function(d,i) {return h-10; })
			.attr("fill", function(d){ return d.color });
	}
}

function addTimeSeries(graph, name, pts, xaxis, yaxis, color) {
	var line = 
		d3.svg.line()
			.x(function(d) { return xaxis( d.x ) } )
			.y(function(d) { return yaxis( d.y ) } )

	graph.selectAll('circle.'+name)
		.data(pts)
	.enter().append('circle')
		.attr("class", "point")
		.attr("r", 3)
		.attr("cx", function(d) { return xaxis( parseFloat(d.x) ) })
		.attr("cy", function(d) { return yaxis( parseFloat(d.y) ) });
	
	graph.selectAll('path.'+name)
		.data([pts])
	.enter().append("path")
		.attr("class", "series")
    	.attr("d", line)
    	.style("stroke", color);
}

function addBarSeries(graph, name, pts, xaxis, yaxis, color) {
	xsize = Array.max(xaxis.range()) - Array.min(xaxis.range());
	barw = xsize / (2 * xaxis.domain().length) - 5;
	
	graph.selectAll('rect.'+name)
			.data(pts)
		.enter().append('rect')
			.attr("class", "bar")
			.attr("x", function(d) { return xaxis(d.x) - barw })
			.attr("y", function(d) { return yaxis(d.y) })
			.attr("width", function(d) { return 2 * barw })
			.attr("height", function(d) { return yaxis(0) - yaxis(d.y) })
			.style("fill", color);
}

function addErrorBars(graph, name, errorbars, xaxis, yaxis, color) {
	ebar = 
		graph.selectAll('g.'+name)
				.data(errorbars)
			.enter().append("g")
				.attr("class", "errorbar")
				.attr("transform", function(d) { return "translate({0},{1})".format(xaxis( d.x ), yaxis( d.y )) });
	
	ebar.append("line")
		.attr("x1", 0)
		.attr("x2", 0)
		.attr("y1", function(d){ return yaxis(d.y + d.s) - yaxis(d.y) })
		.attr("y2", function(d){ return yaxis(d.y - d.s) - yaxis(d.y) })
		.style("stroke", color);
		
	ebar.append("line")
		.attr("x1", -5)
		.attr("x2", 5)
		.attr("y1", function(d){ return yaxis(d.y + d.s) - yaxis(d.y) })
		.attr("y2", function(d){ return yaxis(d.y + d.s) - yaxis(d.y) })
		.style("stroke", color);
	
	ebar.append("line")
		.attr("x1", -5)
		.attr("x2", 5)
		.attr("y1", function(d){ return yaxis(d.y - d.s) - yaxis(d.y) })
		.attr("y2", function(d){ return yaxis(d.y - d.s) - yaxis(d.y) })
		.style("stroke", color);
	
}

function addAxes(graph, xlabels, ylabels, xaxis, yaxis) {
	
	xrange = xaxis.range()
	yrange = yaxis.range()
	
	graph.append('svg:line')
		.attr('class', 'axis')
		.attr('y1', Array.min(yrange) )
		.attr('y2', Array.max(yrange) )
		.attr('x1', Array.min(xrange) )
		.attr('x2', Array.min(xrange) );

	graph.append('svg:line')
		.attr('class', 'axis')
		.attr('y1', Array.max(yrange) )
		.attr('y2', Array.max(yrange) )
		.attr('x1', Array.min(xrange) )
		.attr('x2', Array.max(xrange) );
	

	ticks = 
		graph.selectAll('.ytick')
			.data(ylabels)
		.enter().append('svg:g')
			.attr('transform', function(d) { return "translate("+xaxis(0)+", "+yaxis(d)+")" })
			.attr('class', 'tick');

	ticks.append('svg:line')
		.attr('y1', 0)
		.attr('y2', 0)
		.attr('x1', 0)
		.attr('x2', 5);

	ticks.append('svg:text')
		.text(function(d) { return d; })
		.attr('text-anchor', 'end')
		.attr('dy', 2)
		.attr('dx', -4);
	
	ticks = 
		graph.selectAll('.xtick')
			.data(xlabels)
		.enter().append('svg:g')
			.attr('transform', function(d) { return "translate("+xaxis(d)+", "+yaxis(0)+")" })
			.attr('class', 'tick');
	
	ticks.append('svg:line')
		.attr('y1', -5)
		.attr('y2', 0)
		.attr('x1', 0)
		.attr('x2', 0);
	
	ticks.append('svg:text')
		.text(function(d) { return d; })
		.attr('text-anchor', 'end')
		.attr('dy', 12)
		.attr('dx', 4);
}

function addAxisTitle(graph, title, xaxis, yaxis) {
	xrange = xaxis.range()
	
	xm = ( Array.max(xrange) + Array.min(xrange) ) / 2
	
	graph.append('svg:text')
		.text(title)
		.attr('class', "label")
		.attr('x', xm)
		.attr('y', yaxis(0) + 25 )
		.attr('dy', 6);
}

function createGraph(parent, title, w, h, margin) {
	graph = 
		parent 
			.append("span")
			.append("svg")
			.attr("width", w)
			.attr("height", h)
			.attr("margin", margin);

	graph.append('svg:text')
		.text(title)
		.attr('class', "title")
		.attr('x', w/2)
		.attr('y', margin/2);
	
	return graph;
}