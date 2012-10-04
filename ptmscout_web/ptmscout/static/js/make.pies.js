$(document).ready(function() {
	var num = 1;
	
	d3.selectAll(".data_chart")
		.each(function() {
			json_base64 = 
				d3.select(this)
					.select(".data")
					.text();
			
			data = JSON.parse(Base64.decode(json_base64));
			
			ndata = []
			for(i = 0; i < data.length && i < 10; i++){
				elem = data[i];
				data[i] = {};
				data[i].label = elem[0];
				data[i].x = elem[1];
				ndata.push(data[i])
			}
			
			var div = d3.select(this).append('div');
			createPieChart(div, ndata, 450, 400);
			addExport(div, d3.select(this));
			num+=1;
		});
});