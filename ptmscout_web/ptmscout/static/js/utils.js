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
