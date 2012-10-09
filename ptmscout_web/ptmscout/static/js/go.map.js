function compileNode(name, data){
	val = 0;
	for(i = 0; i < data.length; i++){
		val+=data[i].value;
	}
	return {'GO':name, 'term':"", 'value':val, 'children':data};
}

function createGOMap(json_data, container){
	var w = 1000,
	    h = 800,
	    r = 720,
	    x = d3.scale.linear().range([0, r]),
	    y = d3.scale.linear().range([0, r]),
	    node,
	    root;
	
	var foreground_color = "#000";
	var background_color = "#BBB";
	var foreground_stroke = "none";
	var background_stroke = "none";
	
	var view_limit = 60;
	var id_size = 36;
	var term_size = 24;
	
	var pack = d3.layout.pack()
	    .size([r, r])
	    .value(function(d) { return d.value; })
	
	var vis = container.append("svg:svg")
		.attr("class", "GO")
	    .attr("width", w)
	    .attr("height", h)
	  .append("svg:g")
	    .attr("transform", "translate(" + (w - r) / 2 + "," + (h - r) / 2 + ")");
	

	F_node = compileNode("Molecular Function", json_data.F);
	P_node = compileNode("Biological Process", json_data.P);
	C_node = compileNode("Cellular Component", json_data.C);
	
	root = {'GO':"", 'term':"", 'value':json_data.total, 'children':[F_node, P_node, C_node]};
	
	node = root;

	var nodes = pack.nodes(root);

	vis.selectAll("circle")
    		.data(nodes)
    	.enter().append("svg:circle")
    		.attr("class", function(d) { return d.children ? "parent" : "child"; })
    		.attr("cx", function(d) { return d.x; })
    		.attr("cy", function(d) { return d.y; })
    		.attr("r", function(d) { return d.r; })
    		.on("click", function(d) { return zoom(node == d ? root : d); });
	
	var isNodeVisible = function(d, k) {
		if(k===undefined)
			k = 1;
		return d.children.length != 1 && k * d.r > view_limit && k * d.r <= r / 2;
	};
	
	var isNodeForeground = function(d, parent){
		return d == parent || $.inArray(d, root.children) > -1;
	};
	
	
	var getScaledSize = function(d, size, scale) {
		if(scale == undefined)
			return 2 * d.r / r * size;
		return (2 * (scale(d.y + d.r) - scale(d.y)) / r * size);
	};
	
	vis.selectAll("text.id")
    		.data(nodes.reverse())
    	.enter()
    		.append("svg:text")
	    		.attr("class", function(d) { return "id " + (d.children ? "parent" : "child"); })
	    		.attr("x", function(d) { return d.x; })
	    		.attr("y", function(d) { return d.y; })
	    		.attr("dy", ".35em")
	    		.attr("text-anchor", "middle")
	    		.style("fill", function(d) { return isNodeForeground(d) ? foreground_color : background_color})
	    		.style("stroke", function(d) { return isNodeForeground(d) ? foreground_stroke : background_stroke})
	    		.style("opacity", function(d) { return isNodeVisible(d) ? 1 : 0; })
	    		.style("font-size", function(d) { return getScaledSize(d, id_size) + "px" })
	    		.text(function(d) { return d.GO; });

	vis.selectAll("text.term")
			.data(nodes)
		.enter()
			.append("svg:text")
				.attr("class", function(d) { return "term " + (d.children ? "parent" : "child"); })
	    		.attr("x", function(d) { return d.x; })
	    		.attr("y", function(d) { return d.y; })
	    		.attr("dy", function(d) { return getScaledSize(d, id_size) })
	    		.attr("text-anchor", "middle")
	    		.style("fill", function(d) { return isNodeForeground(d) ? foreground_color : background_color})
	    		.style("stroke", function(d) { return isNodeForeground(d) ? foreground_stroke : background_stroke})
	    		.style("opacity", function(d) { return isNodeVisible(d) ? 1 : 0; })
	    		.style("font-size", function(d) { return getScaledSize(d, term_size) + "px" })
	    		.text(function(d) { return d.term; });

	d3.select("body").on("click", function() { zoom(root); });
	
	
	function zoom(d, i) {
	  var k = r / d.r / 2;
	  x.domain([d.x - d.r, d.x + d.r]);
	  y.domain([d.y - d.r, d.y + d.r]);
	  var parent = d;
	  var t = vis.transition()
	      .duration(d3.event.altKey ? 7500 : 750);
	
	  t.selectAll("circle")
	      .attr("cx", function(d) { return x(d.x); })
	      .attr("cy", function(d) { return y(d.y); })
	      .attr("r", function(d) { return k * d.r; });
	
	  t.selectAll("text.id")
	      .attr("x", function(d) { return x(d.x); })
	      .attr("y", function(d) { return y(d.y); })
	      .style("fill", function(d) { return isNodeForeground(d, parent) ? foreground_color : background_color})
	      .style("stroke", function(d) { return isNodeForeground(d, parent) ? foreground_stroke : background_stroke})
	      .style("opacity", function(d) { return isNodeVisible(d, k) ? 1 : 0; })
	      .style("font-size", function(d) { return getScaledSize(d, id_size, y) + "px" });
	  
	  t.selectAll("text.term")
	      .attr("x", function(d) { return x(d.x); })
	      .attr("y", function(d) { return y(d.y); })
	      .attr("dy", function(d) { return getScaledSize(d, id_size, y); })
	      .style("fill", function(d) { return isNodeForeground(d, parent) ? foreground_color : background_color})
	      .style("stroke", function(d) { return isNodeForeground(d, parent) ? foreground_stroke : background_stroke})
	      .style("opacity", function(d) { return isNodeVisible(d, k) ? 1 : 0; })
	      .style("font-size", function(d) { return getScaledSize(d, term_size, y) + "px" });

	  
	  node = d;
	  d3.event.stopPropagation();
	}
}
	    