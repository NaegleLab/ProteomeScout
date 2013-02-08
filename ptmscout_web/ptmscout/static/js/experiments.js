function toggleVisible(element) {
	experiment_id = element.id;
	$(".subexperiment#"+experiment_id).toggle();
	if($(element).text() === "[+]")
		$(element).text("[-]");
	else
		$(element).text("[+]")
}

function timedRefresh(timeoutPeriod) {
    setTimeout("location.reload(true);",timeoutPeriod);
}

$(document).ready(
	function(){
		$("a.expand").click(
				function(){
					toggleVisible(this);
					return false;
				});
		$("tr.subexperiment").hide();
		
        if( $("div.progress").length > 0 ){
            timedRefresh(10000);
        }
	});
