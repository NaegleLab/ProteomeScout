$(document).ready(function(){
    $(document).tooltip({ track: true });
	json_base64 = d3.select(".GO_map").select(".data").text()
	data = JSON.parse(Base64.decode(json_base64));
	container = d3.select(".GO_map").append("div");
	createGOMap(data, container);
	
	var max_row_display = 10;
	
	addExport(container, d3.select(".GO_map"));
	
	d3.selectAll(".GO_table")
		.each(function(){
			makeTableCollapsable(d3.select(this), 3, max_row_display);
		});
});
