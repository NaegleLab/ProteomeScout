$(document).ready(function(){
	$('.pubinfo').css('display', "none")
	
	$('#new_dataset')
		.click(function(){
			$('#expselect option[value=""]').attr('selected','true')
			$('#parent_exp').css('display', "none");
			$('#change_desc').css('display', "none");
		});
	
	$('#append_dataset')
		.click(function(){
			$('#parent_exp').css('display', "table-row");
			$('#change_desc').css('display', "none");
		});
	
	$('#reload_dataset')
		.click(function(){
			$('#parent_exp').css('display', "table-row");
			$('#change_desc').css('display', "none");
		});
	
	$('#extend_dataset')
		.click(function(){
			$('#parent_exp').css('display', "table-row");
			$('#change_desc').css('display', "table-row");
		});
	
	$('#pubselect')
		.change(function(){
			value = $('#pubselect option:selected').val();
			if(value == "yes")
				$('.pubinfo').css('display',"table-row");
			else
				$('.pubinfo').css('display',"none");
		});
	
	$('#new_dataset').click()
});