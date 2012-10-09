$(document).ready(function(){
	json_base64 = d3.select(".GO_map").select(".data").text()
	data = JSON.parse(Base64.decode(json_base64));
	container = d3.select(".GO_map").append("div");
	createGOMap(data, container);
});