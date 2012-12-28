function update_display() {
    t = window.protein_viewer.transition()
	      .duration(250);

    t.selectAll('rect.residue')
        .each(function(d){
            enabled = 0;
            mod_types = d3.keys(d.value.mods);
            for(var i in mod_types) {
                var k = mod_types[i];
                num_evidences = 0;
                for(var j in d.value.mods[k]){
                    exp_id = d.value.mods[k][j].experiment
                    if(window.experiment_display_modes[exp_id]) {
                        num_evidences+=1;
                    }
                }
                enabled += window.ptm_display_modes[k] && num_evidences > 0 ? 1 : 0;
            }
            d.value.num_mods = enabled;
        })
        .attr('y', function(d) { return window.yaxis(d.value.num_mods); })
        .attr('height', function(d) { return window.yaxis(0) - window.yaxis( d.value.num_mods ); });
}

function toggle_ptm(ptm_name, mode) {
    window.ptm_display_modes[ ptm_name ] = mode;
    update_display();
}

function toggle_exp(exp_id, mode){
    window.experiment_display_modes[exp_id] = mode;
    update_display();
}

function toggle_track(track_name, mode){
    t = window.protein_viewer.transition().duration(250);
    t.selectAll('rect.'+track_name)
        .style('opacity', mode ? 1.0 : 0);
}

function generate_ticks(h, values, cls, show){
    window.protein_viewer.selectAll('line.'+cls)
        .data( values )
            .enter().append('line')
                .attr('class', cls)
                .attr('x1', function(d) { return window.axis(d); })
                .attr('x2', function(d) { return window.axis(d); })
                .attr('y1', window.height/2)
                .attr('y2', window.height/2+h)
                .style('opacity', show ? 1 : 0);

    window.protein_viewer.selectAll('text.'+cls)
        .data( values )
            .enter().append('text')
                .attr('class', cls)
                .attr('x', function(d) { return window.axis(d); })
                .attr('y', window.height/2 + h)
                .attr('dy', '1em')
                .style('opacity', show ? 1 : 0)
                .text(function(d) { return "" + d; });
}

function create_structure_viewer(protein_data) {
    var width = 900;
    var height = 400;
    var domain_height=20;
    var dom_offset = 30;

    window.zoom_level = 1;
    window.height = height;
    window.width = width;


    console.log(protein_data)

    var viewer = 
            d3.select('.protein_viewer .viewer')
                .append('svg')
                    .attr('width', width)
                    .attr('height', height)
                .append('g')
                    .attr('transform', 'translate(0,0)scale(1,1)');

    window.protein_viewer = viewer;
    
    window.experiment_display_modes = {}
    for(var exp_id in protein_data.exps){
        window.experiment_display_modes[exp_id] = true
    }

    window.ptm_display_modes = {}
    for (var i in protein_data.mod_types){
        window.ptm_display_modes[protein_data.mod_types[i]] = true
    }

    var max_mods = 0
    for(var k in protein_data.mods) {
        var num_mods = d3.keys(protein_data.mods[k].mods).length;
        protein_data.mods[k].num_mods = num_mods;
        if(num_mods > max_mods){
            max_mods = num_mods;
        }
    }

    var axis = d3.scale.linear().domain([0, protein_data.seq.length]).range([0, width]);
    var yaxis = d3.scale.linear().domain([0, max_mods]).range([height/2, height/4]);

    window.yaxis = yaxis;
    window.axis = axis;

    var domain_colors = d3.scale.category20();
    var residue_colors = d3.scale.category20();

    viewer
        .append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', width)
            .attr('y1', height/2)
            .attr('y2', height/2);

    viewer
        .append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', width)
            .attr('y1', height/2 + dom_offset)
            .attr('y2', height/2 + dom_offset);

    viewer.selectAll('rect.domain')
        .data(protein_data.domains)
            .enter().append('rect')
                .attr('class', 'domain')
                .attr('x', function(d) { return axis(d.start); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start); })
                .attr('y', height/2 + dom_offset)
                .attr('height', domain_height)
                .attr('title', function(d) { return d.label; })
                .style('fill', function(d) { return domain_colors( d.label ); } )

    viewer.selectAll('rect.residue')
        .data( d3.entries(protein_data.mods) )
            .enter().append('rect')
                .attr('class', 'residue')
                .attr('id', function(d) { return d.value.residue + d.key; })
                .attr('x', function(d) { return axis(d.key - 1); })
                .attr('width', function(d) { return axis(1); })
                .attr('y', function(d) { return yaxis( d.value.num_mods ); })
                .attr('height', function(d) { return yaxis(0) - yaxis( d.value.num_mods ); })
                .style('fill', function(d) { return residue_colors(d.value.residue); });

    ticks = [];
    bigticks = [];
    for(var i = 0; i < protein_data.seq.length; i+=10){
        if(i % 50 == 0)
            bigticks.push(i);
        else
            ticks.push(i);
    }

    generate_ticks(5, ticks, 'tick', protein_data.seq.length < 100);
    generate_ticks(5, bigticks, 'bigtick', true);
}

function get_browser(){
    if (jQuery.browser.mozilla) 
        return '-moz'
    if (jQuery.browser.webkit)
        return '-webkit'
}

$(function(){
    $( document ).tooltip({
                    track: true
                });

    $('.zoomin-tool').button({ icons: { primary: 'ui-icon-zoomin' }, text:false })
                    .click(function(){
                        window.tool = 'z+';
                        $('.viewer').css( 'cursor', get_browser() + '-zoom-in' );
                    });

    $('.zoomout-tool').button({ icons: { primary: 'ui-icon-zoomout' }, text:false })
                    .click(function(){
                         window.tool = 'z-';
                        $('.viewer').css( 'cursor', get_browser() + '-zoom-out' );
                    });

    $('.scroll-tool').button({ icons: { primary: 'ui-icon-arrow-4' }, text:false })
                    .click(function(){
                        window.tool = 's';
                        $('.viewer').css( 'cursor', 'move' );
                    });

    $('.ptm-tool').button()
                  .click(function(){
                    $('.mods').dialog() 
                  });
    $('.exp-tool').button()
                  .click(function(){
                    $('.exps').dialog() 
                  });
    $('.track-tool').button()
                  .click(function(){
                    $('.tracks').dialog() 
                  });
    $('.mods').toggle();
    $('.exps').toggle();
    $('.tracks').toggle();

    $('.tracks #PTMs').change( function() {
        mode = $(this).is(':checked');
        toggle_track('residue', mode);
    });
    $('.tracks #Domains').change( function() {
        mode = $(this).is(':checked');
        toggle_track('domain', mode);
    });

    $('.mods input.modtoggle').change(
        function(){
            mode = $(this).is(':checked');
            ptm = $(this).attr('id');
            toggle_ptm(ptm, mode);
        });

    $('.exps input.exptoggle').change(
        function(){
            mode = $(this).is(':checked');
            exp = $(this).attr('id').substring(1);
            toggle_exp( parseInt(exp), mode );
        });

    $('.protein_viewer').each( function() {
        data = Base64.decode( $(this).find('.data').text() );
        create_structure_viewer( JSON.parse( data ));
    });
});
