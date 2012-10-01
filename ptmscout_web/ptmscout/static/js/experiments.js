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
	});