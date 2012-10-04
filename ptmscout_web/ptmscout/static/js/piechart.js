function createPieChart(node, data, w, h, mouseover) {
	colors = d3.scale.category20();
	
	chart =
		node.append("svg")
			.attr("class", "piechart")
			.attr("width", w)
			.attr("height", h);
	
	origin = 
		chart.append("g")
			.attr("transform", "translate({0},{1})".format(w/2, h/2));
	
	mindim = h;
	if(w < h){
		mindim = w;
	}
	
	dx = 0.0;
	for(i = 0; i < data.length; i++){
		data[i].dx = dx;
		dx += data[i].x
	}
	total = dx;
	
	outlineSize = 3;
	radius = mindim * 0.5 - outlineSize;
	
	angle = function(d, percentage) {
		return 2 * Math.PI * ((d.dx + d.x*percentage) / total);
	};
	
	var arc = d3.svg.arc()
	     .startAngle(function(d) { return angle(d, 0)})
	     .endAngle(function(d) { return angle(d, 1.0)})
	     .innerRadius(0)
	     .outerRadius(radius);
	
	origin.append("circle")
			.attr("class", "outline")
			.attr("cx", 0)
			.attr("cy", 0)
			.attr("r", radius + outlineSize);
	
	origin.selectAll("path.category")
			.data(data)
		.enter().append("path")
			.attr("class", "category")
			.attr("d", arc)
			.style("fill", function(d,i) { return colors(i) });
	
	flip = function(d) {
		if(angle(d, 0.5) > Math.PI && angle(d, 0.5) < Math.PI * 2)
			return -1;
		return 1;
	}
	
	origin.selectAll("g.category_name")
			.data(data)
		.enter()
			.append("g")
				.attr("class", "category_name")
				.attr("transform", function(d) { return "translate({0},{1})scale({3},{4})rotate({2})".format(radius * Math.cos(angle(d, 0.5) - Math.PI/2), radius * Math.sin(angle(d, 0.5) - Math.PI/2), 180 / Math.PI * (angle(d, 0.5) - Math.PI/2), flip(d), flip(d)) })
			.append("text")
				.attr("x", 0)
				.attr("y", 0)
				.attr("dx", function(d){ return - flip(d) * 10 })
				.attr("text-anchor", function(d) { 
					if(flip(d) == -1)
						return "start"
					return "end"
				})
				.attr("dy", "4")
				.text(function(d) { return d.label });
}