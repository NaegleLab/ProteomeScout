QUnit.config.reorder = false;

test( "enable-zoomin", function() {
    ok( window.tool != "z+", "Zoom disabled" );
    $('.zoomin-tool').click()
    if (jQuery.browser.mozilla) {
        ok( $('.viewer').css( 'cursor' ) == '-moz-zoom-in', "Cursor Set" );
    }
    else if (jQuery.browser.webkit) {
        ok( $('.viewer').css( 'cursor' ) == '-webkit-zoom-in', "Cursor Set" );
    }
    ok( window.tool == "z+", "Zoom enabled" );
});

test( "enable-zoomout", function() {
    ok( window.tool != "z-"  , "Zoom disabled" );
    $('.zoomout-tool').click()
    if (jQuery.browser.mozilla) {
        ok( $('.viewer').css( 'cursor' ) == '-moz-zoom-out', "Cursor Set" );
    }
    else if (jQuery.browser.webkit) {
        ok( $('.viewer').css( 'cursor' ) == '-webkit-zoom-out', "Cursor Set" );
    }
    ok( window.tool == "z-", "Zoom enabled" );
});

test( "enable-scroll", function() {
    ok( window.tool != "s", "Scroll disabled" );
    $('.scroll-tool').click()
    ok( $('.viewer').css( 'cursor' ) == 'move', "Cursor Set" );
    ok( window.tool == "s", "scroll enabled" );
});

test( "draw-domains", function() {
    ok( 2 == $('svg rect.domain').length );
});

test( "draw-ptms", function() {
    equal( $('svg rect.residue').length, 27 );
    equal( $('svg rect#Y518').length, 1 );
    equal( $('svg rect#S856').length, 1 );
    equal( $('svg rect#T857').length, 1 );
});

test( "toggle-ptms", function() {
    $('.tracks #PTMs').click();
    stop();

    setTimeout(function() {
        $('svg rect.residue').each(function(){
            equal( $(this).css('opacity'), '0' );
        });
        start();
        $('.tracks #PTMs').click();
    }, 500);
});

test( "toggle-domains", function() {
    $('.tracks #Domains').click();
    stop();

    setTimeout(function() {
        $('svg rect.domain').each(function(){
            equal( $(this).css('opacity'), '0' );
        });
        start();
        $('.tracks #Domains').click();
    }, 500);
});

test( "toggle-ptm-type", function() {
    $('.mods #Phosphoserine').click();
    stop();

    setTimeout(function() {
        $('svg rect.residue').each(function(){
            if($(this).attr('id')[0] == 'S')
                equal( $(this).attr('height'), 0);
        });
        start();
    }, 500);
});

test( "toggle-exp", function() {
    $('.exps #e1304').click();
    stop();

    setTimeout(function() {
        cnt = 0;
        $('svg rect.residue').each(function(){
            cnt += $(this).attr('height') == 0 ? 0 : 1;
        });
        equal( cnt, 5, "Wrong number of visible residues" );
        start();
    }, 500);
});
