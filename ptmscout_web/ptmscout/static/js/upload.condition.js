$(document).ready(function(){
	$("#add_condition").click(function(){
		$(".condition.hidden").first().removeClass("hidden")
	});
	
	$(".remove_condition").click(function(){
		var parent = $(this).parent()
		parent.children('select').val('')
		parent.children('input[type="text"]').val('')
		parent.remove()
	});
});