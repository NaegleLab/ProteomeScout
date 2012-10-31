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
	
	month_vals = ['january','february','march','april','may','june','july','august','september','october','november','december']
	
	$('#load_pubmed')
		.click(function(){
			query_url = $("#pubmed_query_url").text()
			pmid = $('input[name="pmid"]').val()
			pmiderror = $("#pmiderror")
			
			pmiderror.text("Querying...")
			pmiderror.show()
			
			if(/^[0-9]+$/.test(pmid)){
				$.get(query_url + pmid, function(data) {
					  date = data['DA'];
					  year = date.substring(0,4);
					  month = date.substring(4,6)<<0;
					  month_val = month_vals[month - 1];
					  
					  volume = data['VI'];
					  pages = data['PG'].split("-");
					  journal = data['JT'];
					  authors = data['AU'].join();
					  
					  $('input[name="authors"]').val(authors);
					  $('input[name="journal"]').val(journal);
					  $('select[name="publication_month"]').val(month_val);
					  $('input[name="publication_year"]').val(year);
					  $('input[name="volume"]').val(volume);
					  $('input[name="page_start"]').val(pages[0]);
					  $('input[name="page_end"]').val(pages[1]);
					  
					  pmiderror.hide();
					}, 'json')
					.error(function() { 
							pmiderror.text("Pubmed query failed!")
						});
			}else {
				pmiderror.text("PMID must be an integer!")
			}
			
			
		});
});