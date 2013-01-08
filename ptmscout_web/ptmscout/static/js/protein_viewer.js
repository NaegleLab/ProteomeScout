function ZoomWindow(structure_viewer, svg_container, width, height) {
    this.parent_viewer = structure_viewer;
    this.residue = 0;
    this.width = width;
    this.height = height;
    this.max = structure_viewer.protein_data.seq.length;
    this.minwidth = 20;

    console.log(this.max)

    grab_width = 10;

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
                .append('rect')
                    .attr('class', "zoomwindow")
                    .attr('x', axis(zoomer.residue))
                    .attr('width', axis(width))
                    .attr('y', 0)
                    .attr('height', height+50)
                    .call(drag_behavior);

    this.left_expander =
        svg_container
            .append('rect')
                .attr('class', "handle")
                .attr('x', axis(zoomer.residue))
                .attr('width', grab_width)
                .attr('y', 0)
                .attr('height', height+50)
                .call(resize_left_behavior);

    this.right_expander =
        svg_container
            .append('rect')
                .attr('class', "handle")
                .attr('x', axis(zoomer.residue+width)-grab_width)
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

    this.expand_gradient =
        svg_container.selectAll('path.expander')
                .data([this.polygon_vertices])
            .enter().append('path')
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
            .attr('x', nx);
        zoomer.right_expander
            .attr('x', nx + nw - grab_width);

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

ZoomWindow.prototype.hide = function() {
    this.zoom_window.style('opacity', 0);
    this.expand_gradient.style('opacity', 0);
};

ZoomWindow.prototype.hide = function() {
    this.zoom_window.style('opacity', 0.3);
    this.expand_gradient.style('opacity', 1);
};


function TrackViewer(structure_viewer, svg_container, offset, small_ticks) {
    this.domain_height = 20;
    this.dom_offset = 30;
    this.baseline = offset;
    this.barheight = 100;

    var viewer =
            svg_container
                .append('g')
                    .attr('transform', 'translate(0,{0})scale(1,1)'.format(this.baseline));

    this.viewer = viewer;
    this.parent_viewer = structure_viewer;

    var axis = d3.scale.linear().domain(this.parent_viewer.axis.domain()).range(this.parent_viewer.axis.range());
    this.axis = axis;
    var yaxis = structure_viewer.yaxis;
    var residue_colors = structure_viewer.residue_colors;
    var protein_data = structure_viewer.protein_data;

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
                .style('fill', function(d) { return residue_colors(d.value.residue); });

    ticks = [];
    bigticks = [];
    for(var i = 0; i < protein_data.seq.length; i+=10){
        if(i % 50 == 0)
            bigticks.push(i);
        else
            ticks.push(i);
    }

    this.tick_classes = ['tick', 'bigtick'];
    if(small_ticks){
        this.generate_ticks(viewer, 5, ticks, 'tick');
    }
    this.generate_ticks(viewer, 5, bigticks, 'bigtick');
};

TrackViewer.prototype.hide = function() {
    this.viewer.style('opacity', '0');
};

TrackViewer.prototype.show = function() {
    this.viewer.style('opacity', '1');
};

TrackViewer.prototype.update_view = function(axis, yaxis, hide_small_ticks) {
    this.viewer.selectAll('rect.residue')
                .attr('x', function(d) { return axis(d.key - 1); })
                .attr('width', function(d) { return axis(d.key) - axis(d.key - 1); })
                .attr('y', function(d) { return yaxis( d.value.num_mods ); })
                .attr('height', function(d) { return yaxis(0) - yaxis( d.value.num_mods ); });

    var tick_opacity = 1;
    if(hide_small_ticks){
        tick_opacity = 0;
    }

    this.viewer.selectAll('line.tick')
        .style('opacity', tick_opacity);
    this.viewer.selectAll('text.tick')
        .style('opacity', tick_opacity);

    for (var i in this.tick_classes){
        cls = this.tick_classes[i];
        this.viewer.selectAll('line.'+cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); });

        this.viewer.selectAll('text.'+cls)
                .attr('x', function(d) { return axis(d); });
    }

    this.viewer.selectAll('text.aminoacid')
                .attr('x', function(d,i){ return axis(i+0.5); })
                .style('font-size', Math.min(16, axis(1) - axis(0)));

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
    var axis = this.parent_viewer.axis;
    var yaxis = this.parent_viewer.yaxis;
    var domain_colors = this.parent_viewer.domain_colors;
    var track_viewer = this;
    var protein_data = this.parent_viewer.protein_data
    var width = track_viewer.parent_viewer.width;

    this.viewer
        .append('line')
            .attr('class', "strand")
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
};

TrackViewer.prototype.generate_ticks = function(parent_element, h, values, cls){
    var track_viewer = this;
    var axis = this.parent_viewer.axis;

    parent_element.selectAll('line.'+cls)
        .data( values )
            .enter().append('line')
                .attr('class', cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); })
                .attr('y1', 0)
                .attr('y2', h);

    parent_element.selectAll('text.'+cls)
        .data( values )
            .enter().append('text')
                .attr('class', cls)
                .attr('x', function(d) { return axis(d); })
                .attr('y', h)
                .attr('dy', '1em')
                .text(function(d) { return "" + d; } );
};


function StructureViewer(protein_data) {
    this.protein_data = protein_data
    this.width = 900;
    this.height = 550;

    this.macro_viewer_position = 150;
    this.zoom_viewer_position = 400;

    this.svg_container =
            d3.select('.protein_viewer .viewer')
                .append('svg')
                    .attr('width', this.width)
                    .attr('height', this.height);

    this.experiment_display_modes = {}
    for(var exp_id in protein_data.exps){
        this.experiment_display_modes[exp_id] = true;
    }

    this.ptm_display_modes = {}
    for (var i in protein_data.mod_types){
        this.ptm_display_modes[protein_data.mod_types[i]] = true;
    }

    var max_mods = 0
    for(var k in protein_data.mods) {
        var num_mods = d3.keys(protein_data.mods[k].mods).length;
        protein_data.mods[k].num_mods = num_mods;
        if(num_mods > max_mods){
            max_mods = num_mods;
        }
    }

    console.log(protein_data);
    this.barheight=100;

    this.axis = d3.scale.linear().domain([0, protein_data.seq.length]).range([0, this.width]);
    this.yaxis = d3.scale.linear().domain([0, max_mods]).range([0, -this.barheight]);

    this.domain_colors = d3.scale.category20();
    this.residue_colors = d3.scale.category20();

    this.zoom_window = new ZoomWindow(this, this.svg_container, 50, this.macro_viewer_position);

    this.macro_viewer = new TrackViewer(this, this.svg_container, this.macro_viewer_position, protein_data.seq.length < 100);
    this.zoom_viewer = new TrackViewer(this, this.svg_container, this.zoom_viewer_position, true);
    this.zoom_viewer.show_residues();

    this.zoom_window.set_zoom_viewer(this.zoom_viewer);


    this.macro_viewer.show_domains();

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
    t = this.svg_container.transition().duration(250);
    t.selectAll('rect.'+track_name)
        .style('opacity', mode ? 1.0 : 0);
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
        window.structure_viewer = new StructureViewer( JSON.parse( data ));
    });
});
