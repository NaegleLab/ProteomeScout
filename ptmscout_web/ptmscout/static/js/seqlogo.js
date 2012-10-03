function createSeqlogo(node, data, w, h){
	chart = 
		node.append("svg")
			.attr("width", w)
			.attr("height", h+20);
	
	total = data.total;
	total_columns = data.frequencies.length;
	cw = w / total_columns;
	ch = (h - 20) / 12;
	
	colors = d3.scale.category20();
	
	columns = 
		chart.selectAll("g.column")
			.data(data.frequencies)
		.enter().append("g")
			.attr("class", "column")
			.attr("transform", function(d, i) { return "translate({0},0)".format(i * cw) })
	
	entry = 
		columns.selectAll("g.residue")
				.data(function(d){ return d.f; })
			.enter().append("g")
				.attr("class", "residue")
				.attr("transform", function(d, i) { return "translate(0, {0})scale(1.0,{1})".format( (d[2] / total * h) + 10 , d[1] / total * h / ch) });
	
    entry.append("text")
	    .attr("text-anchor", "middle")
	    .attr("x", cw/2)
	    .attr("y", ch - (ch/50))
	    .attr("class", "arial")
	    .style("font-size", 13.5 * (ch / 10) + "px")
	    .style("fill", function(d, i) { return colors(i) })
	    .text(function(d) { return d[0] });
}