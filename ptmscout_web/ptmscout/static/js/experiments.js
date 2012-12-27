function toggleVisible(element) {
	experiment_id = element.id;
	$(".subexperiment#"+experiment_id).toggle();
	if($(element).text() === "[+]")
		$(element).text("[-]");
	else
		$(element).text("[+]")
}

$(document).ready(
	function(){
		$("a.expand").click(
				function(){
					toggleVisible(this);
					return false;
				});
		$("tr.subexperiment").hide();
		
		$("div.progress").each(function(){
				val = $(this).text();
				$(this).text("");
				items = val.split(" / ");
				n = parseInt(items[0]);
				d = parseInt(items[1]);
				
				$(this).progressbar({ value: n, max: d});
			});
	});