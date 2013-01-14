QUnit.config.reorder = false;

test( "draw-domains", function() {
    ok( 4 == $('svg rect.domain').length );
    ok( 4 == $('svg text.domain').length );
    
    visible = [];
    $('svg text.domain').each( function() {
        if( '1' == $(this).css('opacity'))
            visible.push($(this).text());
    });
    
    equal( visible.length, 2 );
    equal( visible[0], "Pkinase_Tyr" );
    equal( visible[1], "Pkinase_Tyr" );
});

test( "draw-ptms", function() {
    equal( $('svg rect.residue').length, 54 );
    equal( $('svg rect#Y518').length, 2 );
    equal( $('svg rect#S856').length, 2 );
    equal( $('svg rect#T857').length, 2 );
});

test( "draw-ticks", function() {
    equal( $('svg line.t10').length, 83 )
    equal( $('svg line.t50').length, 10 )
    equal( $('svg line.t100').length, 16 )
    equal( $('svg line.t500').length, 2 )
    equal( $('svg line.t1000').length, 2 )
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
        equal( cnt, 10, "Wrong number of visible residues" );
        start();
    }, 500);
});

test( "toggle-zoom", function() {
    cur_height = $('svg').attr('height');
    $('.zoomin-tool').click()
    stop();

    setTimeout(function() {
        equal( $('svg rect.window').length, 1 );
        equal( $('svg rect.handle').length, 2 );
        equal( $('svg g.zoom_track_viewer').css('opacity'), 1 );

        new_height = $('svg').css('height');
        ok( cur_height <= new_height );

        start();
    }, 1000);
});

test( "toggle-zoom-off", function() {
    cur_height = $('svg').css('height');
    $('.zoomout-tool').click()
    stop();

    setTimeout(function() {
        equal( $('svg rect.window').length, 0 );
        equal( $('svg rect.handle').length, 0 );
        equal( $('svg g.zoom_track_viewer').css('opacity'), 0 );

        new_height = $('svg').css('height');
        ok( cur_height >= new_height );

        start();
    }, 1000);
});
