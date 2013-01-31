function ZoomWindow(structure_viewer, svg_container, residue, width, height) {
    this.parent_viewer = structure_viewer;
    this.residue = residue;
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

    this.right_expander =
        svg_container
            .insert('rect', ":first-child")
                .attr('class', "zoom handle")
                .attr('x', axis(zoomer.residue+width)-grab_width/2)
                .attr('width', grab_width)
                .attr('y', 0)
                .attr('height', height+50)
                .call(resize_right_behavior);

    this.left_expander =
        svg_container
            .insert('rect', ":first-child")
                .attr('class', "zoom handle")
                .attr('x', axis(zoomer.residue) - grab_width/2)
                .attr('width', grab_width)
                .attr('y', 0)
                .attr('height', height+50)
                .call(resize_left_behavior);


    this.zoom_window =
        svg_container
                .insert('rect', ":first-child")
                    .attr('class', "zoom window")
                    .attr('x', axis(zoomer.residue))
                    .attr('width', axis(width))
                    .attr('y', 0)
                    .attr('height', height+50)
                    .call(drag_behavior);

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
