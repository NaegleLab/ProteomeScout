$(function(){
   $("div.compare_result").each(function(){
        if($(this).hasClass('experiment_result')){
            result_table=$(this).find("table.protein_list");
            $(this).find("button.show_results")
                    .on('click', function(){
                       if(result_table.hasClass('hidden')){
                        result_table.removeClass('hidden');
                        result_table.show();
                       } else {
                        result_table.hide();
                        result_table.addClass('hidden');
                       }
                    });
        }
   }); 
    
});
