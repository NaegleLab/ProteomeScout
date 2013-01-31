QUnit.config.reorder = false;

test( "draw-domains", function() {
    ok( 2 == $('svg rect.domain').length );
    ok( 2 == $('svg text.domain').length );
    
    visible = [];
    $('svg text.domain').each( function() {
        if( '1' == $(this).css('opacity'))
            visible.push($(this).text());
    });
    
    equal( visible.length, 2 );
    equal( visible[0], "Pkinase" );
    equal( visible[1], "Pkinase" );
});

test( "draw-ptms", function() {
    equal( $('svg rect.ptm').length, 72 );
    equal( $('svg rect#T14').length, 2 );
    equal( $('svg rect#Y15').length, 2 );
    equal( $('svg rect#T161').length, 2 );
    ok( -1 != $('#T14').attr('title').indexOf('IEKIGEGtYGVVYKG') );
    ok( -1 != $('#Y15').attr('title').indexOf('EKIGEGTyGVVYKGR') );
    ok( -1 != $('#T161').attr('title').indexOf('GIPIRVYtHEVVTLW') );
});

test( "draw-ticks", function() {
    equal( $('svg line.t10').length, 24 )
    equal( $('svg line.t50').length, 6 )
    equal( $('svg line.t100').length, 4 )
    equal( $('svg line.t500').length, 0 )
    equal( $('svg line.t1000').length, 0 )
});

test( "toggle-ptms", function() {
    $('.tracks #PTMs').click();
    stop();

    setTimeout(function() {
        $('svg g.PTMs').each(function(){
            equal( $(this).css('opacity'), '0' );
        });
        start();
        $('.tracks #PTMs').click();
    }, 1000);
});

test( "toggle-domains", function() {
    $('.tracks #Domains').click();
    stop();

    setTimeout(function() {
        $('svg g.Domains').each(function(){
            equal( $(this).css('opacity'), '0' );
        });
        start();
        $('.tracks #Domains').click();
    }, 1000);
});

test( "toggle-ptm-type", function() {
    $('.mods #Phosphoserine').click();
    stop();

    setTimeout(function() {
        $('svg rect.ptm').each(function(){
            if($(this).attr('id')[0] == 'S')
                equal( $(this).attr('height'), 0);
        });
        start();
    }, 1000);
});

test( "toggle-exp", function() {
    $('.exps #e1304').click();
    stop();

    setTimeout(function() {
        cnt = 0;
        $('svg rect.ptm').each(function(){
            cnt += $(this).attr('height') == 0 ? 0 : 1;
        });
        equal( cnt, 8, "Wrong number of visible residues" );
        start();
    }, 700);
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
