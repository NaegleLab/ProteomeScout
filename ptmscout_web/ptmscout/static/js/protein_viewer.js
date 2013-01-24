function scrollBottom(element){
    var extra = 10;

    var top_pos = element.offset().top;
    var ptm_height = element.height();
    var scrollback = window.innerHeight - ptm_height - extra;

    if( scrollback < 0 )
        scrollback = 0;

    $('html, body').animate({ scrollTop: top_pos - scrollback }, 500);

}

function build_ptm_table(k, mods, protein_data) {
    ms_entries = [];
    for(m in mods.mods) {
        for(i in mods.mods[m]){
            ms = mods.mods[m][i];
            ms.site_pos = k;
            ms.site = "{0}{1}".format(mods.residue, k);
            ms.experiment_name = protein_data.exps[ms.experiment];
            ms.mod_type = m;
            ms.peptide = mods.peptide;
            ms_entries.push(ms);
        }
    }

    d3.select('div.ptm_metadata table').remove();

    table =
        d3.select('div.ptm_metadata').append('table')
                .attr('id', 'peptide_table');

    header_row = table.append('tr');
    header_row.append('th')
            .text('Site');
    header_row.append('th')
            .text('Peptide');
    header_row.append('th')
            .text('Modification');
    header_row.append('th')
            .text('Experiment');
    header_row.append('th')
            .text('Data');

    table.selectAll('tr.moditem')
        .data(ms_entries)
            .enter().append('tr')
                .attr('class', 'moditem')
                .each(function(d){
                    d3.select(this).append('td')
                        .attr('class', 'modsite')
                        .text(d.site);
                    d3.select(this).append('td')
                        .attr('class', 'peptide')
                        .text(d.peptide);
                    d3.select(this).append('td')
                        .text(d.mod_type);
                    exp_td = d3.select(this).append('td');

                    if(d.exported == 1){
                        exp_td.append('a')
                            .attr('href', "{0}/{1}".format(protein_data.experiment_url, d.experiment))
                            .attr('target', '_blank')
                            .text(d.experiment_name);
                    }else{
                        exp_td.text(d.experiment_name);
                    }

                    data_td = d3.select(this).append('td');

                    if(d.has_data) {
                        data_td.append('button')
                            .on('click', function() { window.open("{0}?experiment_id={1}&site_pos={2}".format(protein_data.protein_data_url, d.experiment, d.site_pos), '_blank'); })
                            .text('View');
                    }
                });

    scrollBottom($(".ptm_metadata"))
}

function ZoomWindow(structure_viewer, svg_container, width, height) {
    this.parent_viewer = structure_viewer;
    this.residue = 0;
    this.width = width;
    this.height = height;
    this.max = structure_viewer.protein_data.seq.length;
    this.minwidth = 50;

    grab_width = 8;

    var axis = structure_viewer.axis;
    var zoomer = this;

    drag_behavior =
        d3.behavior.drag()
            .on("drag", drag_window);
    resize_left_behavior =
        d3.behavior.drag()
            .on("drag", resize_left);
    resize_right_behavior =
        d3.behavior.drag()
            .on("drag", resize_right);


    this.zoom_window =
        svg_container
                .insert('rect', ":first-child")
                    .attr('class', "zoom window")
                    .attr('x', axis(zoomer.residue))
                    .attr('width', axis(width))
                    .attr('y', 0)
                    .attr('height', height+50)
                    .call(drag_behavior);

    this.left_expander =
        svg_container
            .insert('rect', ":first-child")
                .attr('class', "zoom handle")
                .attr('x', axis(zoomer.residue) - grab_width/2)
                .attr('width', grab_width)
                .attr('y', 0)
                .attr('height', height+50)
                .call(resize_left_behavior);

    this.right_expander =
        svg_container
            .insert('rect', ":first-child")
                .attr('class', "zoom handle")
                .attr('x', axis(zoomer.residue+width)-grab_width/2)
                .attr('width', grab_width)
                .attr('y', 0)
                .attr('height', height+50)
                .call(resize_right_behavior);

    this.polygon_vertices = 
            [{"x":axis(zoomer.residue),"y":height+50},
             {"x":axis(zoomer.residue+width),"y":height+50},
             {"x":structure_viewer.width,"y":height+150},
             {"x":0,"y":height+150},
            ];

    var line =
		d3.svg.line()
			.x(function(d) { return d.x } )
			.y(function(d) { return d.y } )

    var gradient = svg_container.append("svg:defs")
          .append("svg:linearGradient")
              .attr("id", "gradient")
              .attr("x1", "0%")
              .attr("y1", "0%")
              .attr("x2", "0%")
              .attr("y2", "100%")
              .attr("spreadMethod", "pad");

    gradient.append("svg:stop")
            .attr("offset", "0%")
            .attr("stop-color", "#666")
            .attr("stop-opacity", 0.3);

    gradient.append("svg:stop")
            .attr("offset", "100%")
            .attr("stop-color", "#666")
                .attr("stop-opacity", 0);

    this.gradient = gradient;
    this.expand_gradient =
        svg_container.selectAll('path.expander')
                .data([this.polygon_vertices])
            .enter().insert('path', ":first-child")
                .attr('class', "expander")
                .attr('d', line)
                .style("fill", "url(#gradient)");

    function update_window(){
        nx = axis(zoomer.residue);
        nw = axis(zoomer.width);

        zoomer.polygon_vertices[0].x = nx;
        zoomer.polygon_vertices[1].x = nx+nw;

        zoomer.expand_gradient.attr('d', line(zoomer.polygon_vertices));

        zoomer.zoom_window
            .attr('x', nx)
            .attr('width', nw);
        zoomer.left_expander
            .attr('x', nx - grab_width/2);
        zoomer.right_expander
            .attr('x', nx + nw - grab_width/2);

        if(zoomer.track_viewer != undefined){
            zoomer.track_viewer.view_residues(zoomer.residue, zoomer.width);
        }
    }

    function drag_window() {
        nr = axis.invert(d3.event.x) - zoomer.width/2;
        if(nr < 0)
            nr = 0;
        if(nr + zoomer.width > zoomer.max)
            nr = zoomer.max - zoomer.width;
        zoomer.residue = nr;
        update_window();
    };

    function resize_left() {
        nr = axis.invert(d3.event.x);
        if(nr < 0)
            nr = 0;
        if(nr > zoomer.residue + zoomer.width - zoomer.minwidth)
            nr = zoomer.residue + zoomer.width - zoomer.minwidth;

        zoomer.width += zoomer.residue - nr;
        zoomer.residue = nr;
        update_window();
    }

    function resize_right() {
        nr = axis.invert(d3.event.x);
        if(nr < zoomer.residue + zoomer.minwidth)
            nr = zoomer.residue + zoomer.minwidth;
        if(nr > zoomer.max)
            nr = zoomer.max;

        zoomer.width = nr - zoomer.residue;
        update_window();
    }
}

ZoomWindow.prototype.set_zoom_viewer = function(track_viewer) {
    this.track_viewer = track_viewer;
    this.track_viewer.view_residues(this.residue, this.width);
}

ZoomWindow.prototype.remove = function() {
    this.zoom_window.remove();
    this.gradient.remove();
    this.expand_gradient.remove();
    this.left_expander.remove();
    this.right_expander.remove();
}

ZoomWindow.prototype.hide = function() {
    this.zoom_window.style('opacity', 0);
    this.expand_gradient.style('opacity', 0);
};

ZoomWindow.prototype.hide = function() {
    this.zoom_window.style('opacity', 0.3);
    this.expand_gradient.style('opacity', 1);
};

function get_residue_tooltip_text(data) {
    var site = data.key;
    var amino = data.value.residue;
    var mod_types = d3.keys(data.value.mods);
    var peptide = data.value.peptide;

    var text = "Residue: {0}{1}<br>15-mer: {2}<br>Modifications: {3}<br>".format(amino, site, peptide, mod_types.join(', '));

    return text;
}


function TrackViewer(structure_viewer, svg_container, offset, cls, zoomable) {
    this.domain_height = 20;
    this.dom_offset = 30;
    this.baseline = offset;
    this.barheight = 100;

    var viewer =
            svg_container
                .append('g')
                    .attr('class', cls)
                    .attr('transform', 'translate(0,{0})scale(1,1)'.format(this.baseline));

    this.viewer = viewer;
    this.parent_viewer = structure_viewer;

    var axis = d3.scale.linear().domain(this.parent_viewer.axis.domain()).range(this.parent_viewer.axis.range());
    this.axis = axis;
    var yaxis = structure_viewer.yaxis;
    var residue_colors = structure_viewer.residue_colors;
    var protein_data = structure_viewer.protein_data;
    var track_viewer = this;

    viewer
        .append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', structure_viewer.width)
            .attr('y1', 0)
            .attr('y2', 0);


    viewer.selectAll('rect.residue')
        .data( d3.entries(protein_data.mods) )
            .enter().append('rect')
                .attr('class', 'residue')
                .attr('id', function(d) { return d.value.residue + d.key; })
                .attr('x', function(d) { return axis(d.key - 1); })
                .attr('width', function(d) { return axis(1); })
                .attr('y', function(d) { return yaxis( d.value.num_mods ); })
                .attr('height', function(d) { return yaxis(0) - yaxis( d.value.num_mods ); })
                .attr('title', get_residue_tooltip_text)
                .style('fill', function(d) { return residue_colors(d.value.residue); })
                .style('cursor', 'pointer')
                .on('mouseover', function() { track_viewer.track_mouse_over(this, 'residue', true); })
                .on('mouseout', function() {  track_viewer.track_mouse_over(this, 'residue', false); })
                .on('click', function(d) { build_ptm_table(d.key, d.value, protein_data); });

    this.tick_levels = [5000,1000,500,100,50,10];
    ticks = [];
    for(var t in this.tick_levels)
        ticks.push([]);

    for(var i = 10; i < protein_data.seq.length; i+=10){
        for(var j in this.tick_levels){
            if(i % this.tick_levels[j] == 0){
                ticks[j].push(i);
                break;
            }
        }
    }

    for(var j in this.tick_levels){
        this.generate_ticks(viewer, 5, ticks[j], this.tick_levels[j], zoomable);
    }
};

TrackViewer.prototype.hide = function() {
    this.viewer.style('opacity', '0');
};

TrackViewer.prototype.show = function() {
    this.viewer.style('opacity', '1');
};

TrackViewer.prototype.update_view = function(axis, yaxis) {
    this.viewer.selectAll('rect.residue')
                .attr('x', function(d) { return axis(d.key - 1); })
                .attr('width', function(d) { return axis(d.key) - axis(d.key - 1); })
                .attr('y', function(d) { return yaxis( d.value.num_mods ); })
                .attr('height', function(d) { return yaxis(0) - yaxis( d.value.num_mods ); });

    

    for (var i in this.tick_levels){
        size = this.tick_levels[i];
        cls = 't{0}'.format(size);

        tick_opacity = 0;
        if(this.parent_viewer.width / (axis(size) - axis(0)) < 20){
            tick_opacity = 1;
        }

        this.viewer.selectAll('line.'+cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); })
                .style('opacity', tick_opacity);

        this.viewer.selectAll('text.'+cls)
                .attr('x', function(d) { return axis(d); })
                .style('opacity', tick_opacity);
    }

    this.viewer.selectAll('text.aminoacid')
                .attr('x', function(d,i){ return axis(i+0.5); })
                .style('font-size', Math.min(16, axis(1) - axis(0)));

    this.viewer.selectAll('rect.domain')
                .attr('x', function(d) { return axis(d.start); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start); });

    this.viewer.selectAll('text.domain')
                .each(function(d) { d.show = axis(d.stop) - axis(d.start) > d.label.length * 8 })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .attr('x', function(d) { return (axis(d.start) + axis(d.stop)) / 2; });
}

TrackViewer.prototype.show_residues = function() {
    track_viewer = this;
    this.viewer.selectAll('text.aminoacid')
            .data(this.parent_viewer.protein_data.seq)
                .enter().append('text')
                    .attr('class', 'aminoacid')
                    .attr('x', function(d,i){ return track_viewer.axis(i+0.5); })
                    .attr('y', "-0.1em")
                    .attr('text-anchor', 'middle')
                    .text(function(d) { return d; })
                    .style('font-size', track_viewer.axis(1) - track_viewer.axis(0));
};

TrackViewer.prototype.view_residues = function(residue, width) {
    this.axis.domain([residue, residue+width]);
    yaxis = this.parent_viewer.yaxis;
    this.update_view(this.axis, yaxis, width > 200);
};

TrackViewer.prototype.show_domains = function() {
    var axis = this.axis;
    var yaxis = this.parent_viewer.yaxis;
    var domain_colors = this.parent_viewer.domain_colors;
    var track_viewer = this;
    var protein_data = $.extend(true, [], this.parent_viewer.protein_data);
    var width = track_viewer.parent_viewer.width;

    this.viewer
        .append('line')
            .attr('class', "domain")
            .attr('x1', 0)
            .attr('x2', width)
            .attr('y1', track_viewer.dom_offset)
            .attr('y2', track_viewer.dom_offset);

    this.viewer.selectAll('rect.domain')
        .data(protein_data.domains)
            .enter().append('rect')
                .attr('class', 'domain')
                .attr('x', function(d) { return axis(d.start); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start); })
                .attr('y', track_viewer.dom_offset)
                .attr('height', track_viewer.domain_height)
                .attr('title', function(d) { return d.label; })
                .style('fill', function(d) { return domain_colors( d.label ); } )
                .on('mouseover', function() { track_viewer.track_mouse_over(this, 'domain', true); })
                .on('mouseout', function() { track_viewer.track_mouse_over(this, 'domain', false); })
                .on('click', function(d) { window.open(protein_data.pfam_url + d.label, '_blank'); });

    this.viewer.selectAll('text.domain')
        .data(protein_data.domains)
            .enter().append('text')
                .attr('class', 'domain')
                .attr('x', function(d) { return ( axis(d.start) + axis(d.stop) ) / 2; })
                .attr('y', track_viewer.dom_offset+track_viewer.domain_height)
                .attr('dy', '1em')
                .attr('text-anchor', 'middle')
                .each(function(d) { d.show = ( axis(d.stop) - axis(d.start) ) > d.label.length * 8; })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .text(function(d) { return d.label; });
};

TrackViewer.prototype.track_mouse_over = function(d3element, cls, isover){
    var opacity=0;

    if(this.parent_viewer.track_display_modes[cls]) {
        if(isover)
            opacity=0.8;
        else
            opacity=1.0;
    }

    d3.select(d3element)
        .style('opacity', opacity);
}

TrackViewer.prototype.generate_ticks = function(parent_element, h, values, size, zoomable){
    var track_viewer = this;
    var axis = this.parent_viewer.axis;

    cls = 't{0}'.format(size);
    tick_opacity = 0;
    if(this.parent_viewer.width / (axis(size) - axis(0)) < 20){
        tick_opacity = 1;
    }
    if(!zoomable && tick_opacity == 0){
        return;
    }
    parent_element.selectAll('line.'+cls)
        .data( values )
            .enter().append('line')
                .attr('class', cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); })
                .attr('y1', 0)
                .attr('y2', h)
                .style('opacity', tick_opacity);

    parent_element.selectAll('text.'+cls)
        .data( values )
            .enter().append('text')
                .attr('class', cls)
                .attr('x', function(d) { return axis(d); })
                .attr('y', h)
                .attr('dy', '1em')
                .text(function(d) { return "" + d; } )
                .style('opacity', tick_opacity);
};


function get_15mer(k, seq) {
    k = parseInt(k) - 1;
    var start = k - 7;
    var end = k + 8;

    if(start < 0)
        start = 0;
    if(end > seq.length)
        end = seq.length;

    return seq.substring(start, k) + seq.substring(k,k+1).toLowerCase() + seq.substring(k+1, end);
}


function StructureViewer(protein_data) {
    console.log(protein_data)
    this.show_residues_size_limit = 100;
    this.protein_data = protein_data
    this.width = 900;
    this.height = 250;
    this.default_height=250;
    this.zoom_window_height=550;

    this.macro_viewer_position = 150;
    this.zoom_viewer_position = 400;

    this.svg_container =
            d3.select('.protein_viewer .viewer')
                .append('svg')
                    .attr('width', this.width)
                    .attr('height', this.height);

    this.experiment_display_modes = {};
    for(var exp_id in protein_data.exps){
        this.experiment_display_modes[exp_id] = true;
    }

    this.ptm_display_modes = {};
    for (var i in protein_data.mod_types){
        this.ptm_display_modes[protein_data.mod_types[i]] = true;
    }

    this.track_display_modes = {};
    this.track_display_modes['residue'] = true;
    this.track_display_modes['domain'] = true;

    var max_mods = 0;
    for(var k in protein_data.mods) {
        var num_mods = d3.keys(protein_data.mods[k].mods).length;
        var peptide = get_15mer(k, protein_data.seq)
        protein_data.mods[k].num_mods = num_mods;
        protein_data.mods[k].peptide = peptide;
        if(num_mods > max_mods){
            max_mods = num_mods;
        }
    }

    this.barheight=100;

    this.axis = d3.scale.linear().domain([0, protein_data.seq.length]).range([0, this.width]);
    this.yaxis = d3.scale.linear().domain([0, max_mods]).range([0, -this.barheight]);

    this.domain_colors = d3.scale.category20();
    this.residue_colors = d3.scale.category20();

    this.macro_viewer = new TrackViewer(this, this.svg_container, this.macro_viewer_position, 'macro_track_viewer', false);
    this.macro_viewer.show_domains();
    if(protein_data.seq.length < this.show_residues_size_limit){
        this.macro_viewer.show_residues();
    }

    this.zoom_viewer = new TrackViewer(this, this.svg_container, this.zoom_viewer_position, 'zoom_track_viewer', true);
    this.zoom_viewer.show_residues();
    this.zoom_viewer.show_domains();
    this.zoom_viewer.hide();

    this.zoom_enabled = false;
};

StructureViewer.prototype.zoom_off = function(){
    if(this.zoom_enabled){
        this.zoom_enabled=false;

        this.height = this.default_height;
        var viewer = this;

        $('.protein_viewer svg').animate({
                height: viewer.height
          }, 750, function() {
            viewer.zoom_window.remove();
            viewer.zoom_viewer.hide();
          });
    }
}

StructureViewer.prototype.zoom_on = function() {
    if(!this.zoom_enabled){
        this.zoom_enabled=true;

        this.zoom_window = new ZoomWindow(this, this.svg_container, 50, this.macro_viewer_position);
        this.zoom_window.set_zoom_viewer(this.zoom_viewer);
        this.zoom_viewer.show();

        this.height = this.zoom_window_height;

        $('.protein_viewer svg').animate({
                height: this.height
          }, 750, function() {
          });
    }
};

StructureViewer.prototype.update_display = function() {
    var protein_viewer = this;
    t = this.svg_container.transition()
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
                    if(protein_viewer.experiment_display_modes[exp_id]) {
                        num_evidences+=1;
                    }
                }
                enabled += protein_viewer.ptm_display_modes[k] && num_evidences > 0 ? 1 : 0;
            }
            d.value.num_mods = enabled;
        })
        .attr('y', function(d) { return protein_viewer.yaxis(d.value.num_mods); })
        .attr('height', function(d) { return protein_viewer.yaxis(0) - protein_viewer.yaxis( d.value.num_mods ); });


    domain_mode = this.track_display_modes['domain']

    t.selectAll("line.domain")
        .style('opacity', domain_mode ? 1 : 0);
    t.selectAll("rect.domain")
        .style('opacity', domain_mode ? 1 : 0)
        .style('cursor', domain_mode ? 'pointer' : 'default');
    t.selectAll("text.domain")
        .style('opacity', function(d) { return domain_mode && d.show ? 1 : 0; });

    residue_mode = this.track_display_modes['residue'];
    t.selectAll(".residue")
        .style('opacity', residue_mode ? 1 : 0)
        .style('cursor', residue_mode ? 'pointer' : 'default');
}

StructureViewer.prototype.toggle_ptm = function(ptm_name, mode) {
    this.ptm_display_modes[ ptm_name ] = mode;
    this.update_display();
}

StructureViewer.prototype.toggle_exp = function(exp_id, mode){
    this.experiment_display_modes[exp_id] = mode;
    this.update_display();
}

StructureViewer.prototype.toggle_track = function(track_name, mode){
    this.track_display_modes[ track_name ] = mode;
    this.update_display();
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
                         window.structure_viewer.zoom_on();
                      });

    $('.zoomout-tool').button({ icons: { primary: 'ui-icon-zoomout' }, text:false })
                      .click(function(){
                         window.structure_viewer.zoom_off();
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
        window.structure_viewer.toggle_track('residue', mode);
    });
    $('.tracks #Domains').change( function() {
        mode = $(this).is(':checked');
        window.structure_viewer.toggle_track('domain', mode);
    });

    $('.mods input.modtoggle').change(
        function(){
            mode = $(this).is(':checked');
            ptm = $(this).attr('id');
            window.structure_viewer.toggle_ptm(ptm, mode);
        });

    $('.exps input.exptoggle').change(
        function(){
            mode = $(this).is(':checked');
            exp = $(this).attr('id').substring(1);
            window.structure_viewer.toggle_exp( parseInt(exp), mode );
        });

    $('.protein_viewer').each( function() {
        data = Base64.decode( $(this).find('.data').text() );
        json_data = JSON.parse( data );
        window.structure_viewer = new StructureViewer( json_data );
        if(json_data.experiment != null && json_data.experiment != undefined){
            $('.exps input.exptoggle').each(
                function(){
                    exp = parseInt($(this).attr('id').substring(1));
                    if(exp != json_data.experiment)
                        $(this).click();
                });
        }
    });
});
