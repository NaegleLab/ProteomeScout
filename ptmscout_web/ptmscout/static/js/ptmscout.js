$(function(){
    $(".longtask")
        .on('click', function(e) {
            $("#waiting-dialog")
                .dialog({
                    height: 175,
                    width: 150,
                    modal: true,
                    draggable: false,
                    resizable: false
                });
        });

    $("div.progress").each(function(){
            var val = $(this).text();
            $(this).text("");
            var items = val.split(" / ");
            var n = parseInt(items[0]);
            var d = parseInt(items[1]);
            
            if(n == 0){
                $(this).progressbar({ value: false });
            }else{
                $(this).progressbar({ value: n, max: d});
            }
        });


});
