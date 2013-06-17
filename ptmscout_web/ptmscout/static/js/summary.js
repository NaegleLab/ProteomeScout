$(document).ready( function (){
	json_base64 = Base64.decode($(".seqdata").text());
	data = JSON.parse(json_base64);
	createSeqlogo(d3.select(".seqchart"), data, 400, 325);
});
