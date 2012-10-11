function setValues(exp){
	if(exp === undefined){
		exp = {'name':"",'contact':"",'author':"",'journal':"",'volume':"",'page_start':"",'page_end':"", 'publication_month':"", 'publication_year':""};
		published = ""
	} else
		published = (exp.published == 1 ? "yes" : "no");
	
	$('input[name="experiment_name"]').val(exp.name)
	$('input[name="author_contact"]').val(exp.contact)
	$('input[name="pmid"]').val(exp.PMID)
	$('input[name="URL"]').val(exp.URL)
	
	$('select[name="published"] option[value="'+published+'"]').attr('selected','true')
	
	if(published == "yes"){
		$('.pubinfo').css('display',"table-row");
		
		$('input[name="authors"]').val(exp.author)
		$('input[name="journal"]').val(exp.journal)
		$('input[name="volume"]').val(exp.volume)
		$('input[name="page_start"]').val(exp.page_start)
		$('input[name="page_end"]').val(exp.page_end)
		$('select[name="publication_month"] option[value="'+exp.publication_month+'"]').attr('selected','true')
		$('input[name="publication_year"]').val(exp.publication_year)
	} else
		$('.pubinfo').css('display',"none");
	
	
}


function getExperiment(exp_id, experiments){
	for(var i=0; i < experiments.length; i++){
		if(experiments[i].id == exp_id)
			return experiments[i]
	}
}

$(document).ready(function(){
	
	json_base64 = $('#user_experiments').text();
	window.user_data = JSON.parse(Base64.decode(json_base64));
	
	$('#parent_exp').css('display', "none");
	$('#change_desc').css('display', "none");
	$('.pubinfo').css('display', "none")


	$('#new_dataset')
		.click(function(){
			$('#parent_exp').css('display', "none");
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
		
	$('#expselect')
		.change(function(){
			i = $('#expselect option:selected').val();
			if(i != ''){
				exp = getExperiment(parseInt(i), window.user_data)
				setValues(exp);
			}
			else
				setValues();
		});
	
	$('#pubselect')
		.change(function(){
			value = $('#pubselect option:selected').val();
			if(value == "yes")
				$('.pubinfo').css('display',"table-row");
			else
				$('.pubinfo').css('display',"none");
		});
});