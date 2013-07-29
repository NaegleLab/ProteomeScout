function displaySVG(id) {
	var svgml = $(id).html();
    svgml = svgml.replace("#zoomgrad", "url(#zoomgrad)");
    var b64 = Base64.encode(svgml);
    window.open("data:image/svg+xml;base64,\n"+b64);
}

function exportSVG(chart){
	var id = "#" + chart.attr('id');
	
	if($(id + " svg style").size() == 0){
		
		css_url = $("#graph_css_export_url").text()
		
		$.get(css_url,
			function(data) {
				$(id + " svg").append("<style>" + data + "</style>");
				displaySVG(id);
			}
		).error(function() {
			alert("Request ERROR: unable to load stylesheet. Please try again.")
		});
		
	} else {
		displaySVG(id);
	}
}

var export_num = 0;

function addExport(parent, container) {
	var svg = 
		parent.select('svg')
			.attr('version', "1.1")
			.attr('xmlns', "http://www.w3.org/2000/svg");
	
	parent.attr('id', "chart" + export_num);
	export_num++;
	
	container
		.style('text-align', "right")
		.append('div')
			.attr('class', 'absolute-lr')
            .style('z-index', "999")
			.append('button')
				.text("Export SVG")
				.on('click', function() { exportSVG(parent) });
}
