function compileNode(name, data){
	val = 0;
	for(i = 0; i < data.length; i++){
		val+=data[i].value;
	}
	return {'GO':name, 'term':name, 'value':val, 'children':data};
}

function createGOMap(json_data, container){
	var w = 1000,
	    h = 800,
	    r = 720,
	    x = d3.scale.linear().range([0, r]),
	    y = d3.scale.linear().range([0, r]),
	    node,
	    root;
	
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
	
	root = {'GO':"", 'term':"GO annotations", 'value':json_data.total, 'children':[F_node, P_node, C_node]};
	
	console.log(root);
	
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

	vis.selectAll("text")
    		.data(nodes)
    	.enter().append("svg:text")
    		.attr("class", function(d) { return d.children ? "parent" : "child"; })
    		.attr("x", function(d) { return d.x; })
    		.attr("y", function(d) { return d.y; })
    		.attr("dy", ".35em")
    		.attr("text-anchor", "middle")
    		.style("opacity", function(d) { return d.children.length != 1 && d.r > 20 ? 1 : 0; })
    		.style("font-size", function(d) { return d.r / 10 })
    		.text(function(d) { return d.GO; });

	d3.select("body").on("click", function() { zoom(root); });
	
	
	function zoom(d, i) {
	  var k = r / d.r / 2;
	  x.domain([d.x - d.r, d.x + d.r]);
	  y.domain([d.y - d.r, d.y + d.r]);
	
	  var t = vis.transition()
	      .duration(d3.event.altKey ? 7500 : 750);
	
	  t.selectAll("circle")
	      .attr("cx", function(d) { return x(d.x); })
	      .attr("cy", function(d) { return y(d.y); })
	      .attr("r", function(d) { return k * d.r; });
	
	  t.selectAll("text")
	      .attr("x", function(d) { return x(d.x); })
	      .attr("y", function(d) { return y(d.y); })
	      .style("opacity", function(d) { return d.children.length != 1 && k * d.r > 20 && k * d.r <= r / 2 ? 1 : 0; })
	      .style("font-size", function(d) { return (y(d.y + d.r) - y(d.y - d.r)) / 20 });
	
	  node = d;
	  d3.event.stopPropagation();
	}
}
	    