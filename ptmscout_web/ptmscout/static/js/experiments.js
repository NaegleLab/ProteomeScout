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
		
		$("div.progress").each(function(){
				val = $(this).text();
				$(this).text("");
				items = val.split(" / ");
				n = parseInt(items[0]);
				d = parseInt(items[1]);

                if(n == 0){
                    n = false;
                }
				$(this).progressbar({ value: n, max: d});
			});

        if( $("div.progress").length > 0 ){
            timedRefresh(10000);
        }
	});
