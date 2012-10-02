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
				.attr("transform", function(d, i) { return "translate(0, {0})scale(2.0,{1})".format( (d[2] / total * h) + 10 , d[1] / total * h / 15) });
		
/*	
	entry.append("rect")
			.attr("x", 0)
			.attr("y", 0)
			.attr("width", column_width / 2.0)
			.attr("height", function(d) { return 15 })
			.style("fill", function(d, i) { return colors(i) })
	*/
	entry.append("text")
			.attr("x", 0)
			.attr("y", 15)
			.style("font-size", "23px")
			.style("fill", function(d, i) { return colors(i) })
			.text(function(d){ return d[0] });

	
}