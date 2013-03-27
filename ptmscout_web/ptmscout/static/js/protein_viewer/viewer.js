function TrackViewer(structure_viewer, svg_container, offset, cls) {
    this.baseline = offset;
    this.structure_viewer = structure_viewer;
    this.viewer =
            svg_container
                .append('g')
                    .attr('class', cls)
                    .attr('transform', 'translate(0,{0})scale(1,1)'.format(this.baseline));
    this.axis = d3.scale.linear().domain(structure_viewer.axis.domain()).range(structure_viewer.axis.range());

    this.tracks = [];
};

TrackViewer.prototype.hide = function() {
    this.viewer.style('opacity', '0');
};

TrackViewer.prototype.show = function() {
    this.viewer.style('opacity', '1');
};

TrackViewer.prototype.view_residues = function(residue, width) {
    this.axis.domain([residue, residue+width]);

    for(var i in this.tracks) {
        this.tracks[i].update_display(this.axis, this.structure_viewer.width);
    }
};

TrackViewer.prototype.add_track = function(track) {
    this.tracks.push(track);
    pos = this.get_track_position(track.name);

    track.g.attr('transform', 'translate(0,{0})'.format(pos));
}

TrackViewer.prototype.get_track_position = function(track_name) {
    var npos = 0;
    for(var i in this.tracks){
        if(this.tracks[i].name == track_name){
            return npos;
        }else if(this.tracks[i].visible){
            npos += this.tracks[i].height;
        }
    }
}

TrackViewer.prototype.get_track = function(track_name) {
    for(var i in this.tracks){
        if(this.tracks[i].name == track_name){
            return this.tracks[i];
        }
    }
}

TrackViewer.prototype.toggle_track = function(t, track_name, mode) {
    track = this.get_track(track_name);
    var track_pos = this.get_track_position(track_name);
    var ntransform;
    var ctransform;
    var npos = 0;
    var offscreen_pos = (this.structure_viewer.protein_data.seq.length / 50) * this.structure_viewer.width;

    track.visible = mode;

    if(mode){
        track_pos = undefined;

        for(var i in this.tracks){
            var hpos = 0;

            if(this.tracks[i].name == track_name){
                hpos = offscreen_pos;
                track_pos = npos;
            }

            if(this.tracks[i].visible){
                ntransform = 'translate({0}, {1})'.format(hpos, npos);
                t.select(this.tracks[i].selector)
                    .attrTween('transform', d3.tween(ntransform, d3.interpolateTransform));
                npos += this.tracks[i].height;
            }
        }

        t = t.transition().delay(250);
        ntransform = 'translate(0, {0})'.format(track_pos);
        t.select(track.selector)
                    .attrTween('transform', d3.tween(ntransform, d3.interpolateTransform));

        t = t.transition().delay(500);
        t.select(track.selector)
            .style('opacity', 1);

    }else{

        ctransform = track.g.attr('transform');
        ntransform = 'translate({0},{1})'.format(offscreen_pos, track_pos);

        t.select(track.selector)
                .style('opacity', 0);

        t = t.transition().delay(250);

        t.select(track.selector)
            .attrTween('transform', d3.tween(ntransform, d3.interpolateTransform));

        t = t.transition().delay(500);
        npos = 0;
        for(var i in this.tracks){
            if(this.tracks[i].visible){
                ntransform = 'translate(0, {0})'.format(npos);
                t.select(this.tracks[i].selector)
                    .attrTween('transform', d3.tween(ntransform, d3.interpolateTransform));
                npos += this.tracks[i].height;
            }
        }
    }
}


function StructureViewer(protein_data) {
    console.log(protein_data)
    this.protein_data = protein_data;

    this.show_residues_size_limit = 100;
    this.default_height = 340;
    this.macro_viewer_position = 0;
    this.zoom_viewer_position = 340;
    this.zoom_window_height = 640;

    this.transition_duration = 250;

    this.last_zoom_residue = 0;
    this.last_zoom_width = 50;

    this.width = 900;
    this.height = this.default_height;

    this.svg_container =
            d3.select('.protein_viewer .viewer')
                .append('svg')
                    .attr('width', this.width)
                    .attr('height', this.height);


    this.axis = d3.scale.linear().domain([0, protein_data.seq.length]).range([0, this.width]);

    this.domain_colors = d3.scale.category20();
    this.region_colors = d3.scale.category20b();
    this.residue_colors = create_amino_acid_colors();

    var macro_residues = this.protein_data.seq.length <= this.show_residues_size_limit;
    this.macro_viewer = new TrackViewer(this, this.svg_container, this.macro_viewer_position, 'macro_track_viewer', macro_residues);
    this.create_empty_track(this.macro_viewer);
    this.create_ptm_track(this.macro_viewer);
    this.create_residue_track(this.macro_viewer, this.show_residues_size_limit >= this.protein_data.seq.length);

    this.create_mutation_track(this.macro_viewer, macro_residues);

    this.create_region_track(this.macro_viewer);
    this.create_domain_track(this.macro_viewer);

    this.zoom_viewer = new TrackViewer(this, this.svg_container, this.zoom_viewer_position, 'zoom_track_viewer', true);
    this.zoom_viewer.view_residues(this.last_zoom_residue, this.last_zoom_width);

    this.create_empty_track(this.zoom_viewer);
    this.create_ptm_track(this.zoom_viewer);
    this.create_residue_track(this.zoom_viewer, true);

    this.create_mutation_track(this.zoom_viewer, true);

    this.create_region_track(this.zoom_viewer);
    this.create_domain_track(this.zoom_viewer);
    this.zoom_viewer.hide();

    this.zoom_enabled = false;
};

StructureViewer.prototype.create_ptm_track = function(track_viewer) {
    ptm_track = new PTMTrack('PTMs', track_viewer.viewer, this.protein_data);
    ptm_track.create(track_viewer.axis, this.width, this.residue_colors);
    track_viewer.add_track(ptm_track);
};

StructureViewer.prototype.create_mutation_track = function(track_viewer, show_residues) {
    mutation_track = new MutationTrack('Mutations', track_viewer.viewer, this.protein_data);
    mutation_track.create(track_viewer.axis, this.width, show_residues);
    track_viewer.add_track(mutation_track);
};

StructureViewer.prototype.create_residue_track = function(track_viewer, show_residues) {
    residue_track = new ResidueTrack('Residues', track_viewer.viewer, this.protein_data);
    residue_track.create(track_viewer.axis, this.width, show_residues);
    track_viewer.add_track(residue_track);
};

StructureViewer.prototype.create_region_track = function(track_viewer) {
    region_track = new RegionTrack('Regions', track_viewer.viewer, this.protein_data);
    region_track.create(track_viewer.axis, this.width, this.region_colors);
    track_viewer.add_track(region_track);
};

StructureViewer.prototype.create_domain_track = function(track_viewer) {
    domain_track = new DomainTrack('Domains', track_viewer.viewer, this.protein_data);
    domain_track.create(track_viewer.axis, this.width, this.domain_colors);
    track_viewer.add_track(domain_track);
};

StructureViewer.prototype.create_empty_track = function(track_viewer) {
    empty_track = new EmptyTrack('None', track_viewer.viewer, this.protein_data);
    track_viewer.add_track(empty_track);
};

StructureViewer.prototype.zoom_off = function(){
    if(this.zoom_enabled){
        this.zoom_enabled=false;

        this.last_zoom_residue = this.zoom_window.residue;
        this.last_zoom_width = this.zoom_window.width;

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

        this.zoom_window = new ZoomWindow(this, this.svg_container, this.last_zoom_residue, this.last_zoom_width, this.zoom_viewer_position - 100);
        this.zoom_window.set_zoom_viewer(this.zoom_viewer);
        this.zoom_viewer.show();

        this.height = this.zoom_window_height;

        $('.protein_viewer svg').animate({
                height: this.height
          }, 750, function() {
          });
    }
};

StructureViewer.prototype.toggle_ptm = function(ptm_name, mode) {
    track = this.macro_viewer.get_track('PTMs');
    track.toggle_ptm(ptm_name, mode);

    track = this.zoom_viewer.get_track('PTMs');
    track.toggle_ptm(ptm_name, mode);

    this.macro_viewer.get_track('PTMs').update_values( this.transition_duration );
    this.zoom_viewer.get_track('PTMs').update_values( this.transition_duration );
}

StructureViewer.prototype.toggle_exp = function(exp_id, mode){
    track = this.macro_viewer.get_track('PTMs');
    track.toggle_exp(exp_id, mode);

    track = this.zoom_viewer.get_track('PTMs');
    track.toggle_exp(exp_id, mode);

    this.macro_viewer.get_track('PTMs').update_values( this.transition_duration );
    this.zoom_viewer.get_track('PTMs').update_values( this.transition_duration );
}

StructureViewer.prototype.toggle_track = function(track_name, mode){
    t = this.svg_container.transition()
              .duration(this.transition_duration);

    this.macro_viewer.toggle_track(t, track_name, mode);
    this.zoom_viewer.toggle_track(t, track_name, mode);
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

    
    $('.tracks input').change( function() {
        track = $(this).attr('id');
        mode = $(this).is(':checked');
        window.structure_viewer.toggle_track(track, mode);
    });

    $('.mods input.modtoggle').change(
        function(){
            mode = $(this).is(':checked');
            ptm = $(this).attr('id').replace("_"," ");
            window.structure_viewer.toggle_ptm(ptm, mode);
        });

    $('.exps input.exptoggle').change(
        function(){
            mode = $(this).is(':checked');
            exp = $(this).attr('id').substring(1);
            window.structure_viewer.toggle_exp( parseInt(exp), mode );
        });

    $('.mods button.all').click(
        function(){
            $(".mods input.modtoggle").each(
                function(){
                    if(! $(this).is(':checked')) $(this).click();
                });
        });

    $('.mods button.none').click(
        function(){
            $(".mods input.modtoggle").each(
                function(){
                    if($(this).is(':checked')) $(this).click();
                });
        });

    $('.exps button.all').click(
        function(){
            $(".exps input.exptoggle").each(
                function(){
                    if(! $(this).is(':checked')) $(this).click();
                });
        });

    $('.exps button.none').click(
        function(){
            $(".exps input.exptoggle").each(
                function(){
                    if($(this).is(':checked')) $(this).click();
                });
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
