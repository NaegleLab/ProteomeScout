function createSeqlogo(node, data, w, h){
	chart = 
		node.append("svg")
			.attr("width", w)
			.attr("height", h+20);
	
	total = data.total;
	total_columns = data.frequencies.length;
	column_width = w / total_columns;
	
	colors = d3.scale.category20();
	
	columns = 
		chart.selectAll("g.column")
			.data(data.frequencies)
		.enter().append("g")
			.attr("class", "column")
			.attr("transform", function(d, i) { return "translate({0},0)".format(i * column_width) })
	
	entry = 
		columns.selectAll("g.residue")
				.data(function(d){ return d.f; })
			.enter().append("g")
				.attr("class", "residue")
				.attr("transform", function(d, i) { return "translate(0, {0})scale(2.0,{1})".format( (d[2] / total * h) + 10 , 1.04 * d[1] / total * h / 20) });
	
	entry.append("text")
			.attr("x", 0)
			.attr("y", 19)
			.style("font-size", "30px")
			.style("fill", function(d, i) { return colors(i) })
			.text(function(d){ return d[0] });
}