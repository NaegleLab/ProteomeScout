Array.max = function( array ){
    return Math.max.apply( Math, array );
};
Array.min = function( array ){
    return Math.min.apply( Math, array );
};

Array.unique = function( array ){
	return array.filter(function(itm,i,a){
			return i==a.indexOf(itm);
		});
};

String.prototype.format = function() {
	  var args = arguments;
	  return this.replace(/{(\d+)}/g, function(match, number) { 
	    return typeof args[number] != 'undefined'
	      ? args[number]
	      : match
	    ;
	  });
	};


function done_waiting(){
    $(".waiting-modal").hide();
}

$(function(){
    $(".longtask")
        .on('click', function(e) {
            $(".waiting-modal").show();
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
