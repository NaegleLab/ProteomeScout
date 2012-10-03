function exportSVG(chart){
	var id = "#" + chart.attr('id');
	
	var svgml = $(id).html();
    var b64 = Base64.encode(svgml);
    window.open("data:image/svg+xml;base64,\n"+b64)
}

var export_num = 0;

function addExport(parent, container) {
	var svg = 
		parent
			.select('svg')
			.attr('version', "1.1")
			.attr('xmlns', "http://www.w3.org/2000/svg");
	
	parent.attr('id', "chart" + export_num)
	export_num++;
	
	container
		.append('div')
		.style("text-align", "right")
			.append('button')
			.text("Export SVG")
			.on('click', function() { exportSVG(parent) });
}