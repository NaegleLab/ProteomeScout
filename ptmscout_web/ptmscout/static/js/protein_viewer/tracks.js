// Track interface:
//   __init__(name, track_viewer, protein_data)
//   create(axis, width, ...)                      // should create initial elements
//   update_display(axis, width)                   // should update visible elements, zoom, etc
//   update_values(transition_duration)     // should animate value transitions
//   height                                 // the vertical track size for spacing
//   g                                      // the container element

// Bookeeping variables for TrackViewer:
//   name                                   // the name of this track
//   visible                                // the visibility status of this track
//   pos                                    // the current vertical offset position
//   id                                     // the dom node id of the track g container element
//   selector                               // the dom node selector

// Init code:
//   this.track_viewer = track_viewer;
//   this.protein_data = protein_data;
//   this.name = name;
//   this.visible = true;
//   this.g = track_viewer.append('g');


// Track types:
//    Residue track
//    PTM track:
//      toggle_ptm(name, mode)
//      toggle_exp(eid, mode)
//    Domain track
//    Region track
//    Mutation track

window.next_track_id = 0;

function init_track(track, name, track_viewer, protein_data) {
    track.track_viewer = track_viewer;
    track.protein_data = protein_data;
    track.name = name;
    track.visible = true;

    track.id = 'track{0}'.format(window.next_track_id++);
    track.selector = '#{0}'.format(track.id);
    track.g = track_viewer.append('g')
                            .attr('id', track.id)
                            .attr('class', name);

};

function PTMTrack(name, track_viewer, protein_data) {
    init_track(this, name, track_viewer, protein_data);

    // configurables
    this.height = 100;
    this.barheight = 100;
    this.min_barwidth = 1;

    var max_mods = 0;
    for(var k in protein_data.mods) {
        var num_mods = d3.keys(protein_data.mods[k].mods).length;
        protein_data.mods[k].num_mods = num_mods;
        if(num_mods > max_mods){
            max_mods = num_mods;
        }
    }

    this.yaxis = d3.scale.linear().domain([0, max_mods]).range([this.height, this.height-this.barheight]);

    this.g.append('text')
        .attr('x', 0)
        .attr('y', 0)
        .attr('text-anchor', 'left')
        .attr('class', 'track-label')
        .text(name);
};

function scrollBottom(element){
    var extra = 10;

    var top_pos = element.offset().top;
    var ptm_height = element.height();
    var scrollback = window.innerHeight - ptm_height - extra;

    if( scrollback < 0 )
        scrollback = 0;

    $('html, body').animate({ scrollTop: top_pos - scrollback }, 500);

};

function is_kinase_loop(k, protein_data){
    for(var i in protein_data.regions){
        if(protein_data.regions[i].start <= k && k <= protein_data.regions[i].stop){
            if(protein_data.regions[i].label == 'Kinase Activation Loop')
                return 'K';
            if(protein_data.regions[i].label == 'Possible Kinase Activation Loop')
                return '?';
        }
    }
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
            ms.is_mutated = (k in protein_data.mutations);
            ms.kinase_loop = is_kinase_loop(k, protein_data);

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
            .text('Annotations');
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

                    annotation_td = d3.select(this).append('td');
                    images_url = 'http://ptmscout.wustl.edu/static/images/{0}';
                    if(d.is_mutated)
                        annotation_td.append('img')
                                .attr('title', 'This residue has recorded natural variants')
                                .attr('src', images_url.format('red_flag.jpg'));
                    if(d.kinase_loop == 'K')
                        annotation_td.append('img')
                                .attr('title', 'This residue is in a predicted kinase activation loop')
                                .attr('src', images_url.format('kinase.jpg'));
                    if(d.kinase_loop == '?')
                        annotation_td.append('img')
                                .attr('title', 'This residue is in a possible kinase activation loop')
                                .attr('src', images_url.format('question.jpg'));

                    data_td = d3.select(this).append('td');


                    if(d.has_data) {
                        data_td.append('button')
                            .on('click', function() { window.open("{0}?experiment_id={1}&site_pos={2}".format(protein_data.protein_data_url, d.experiment, d.site_pos), '_blank'); })
                            .text('View');
                    }
                });

    scrollBottom($(".ptm_metadata"));
};

PTMTrack.prototype.create = function(axis, viewer_width, residue_colors) {
    ptm_viewer = this;

    this.experiment_display_modes = {};
    for(var exp_id in this.protein_data.exps){
        this.experiment_display_modes[exp_id] = true;
    }

    this.ptm_display_modes = {};
    for (var i in this.protein_data.mod_types){
        this.ptm_display_modes[this.protein_data.mod_types[i]] = true;
    }

    this.g.selectAll('rect.ptm')
        .data( d3.entries(this.protein_data.mods) )
            .enter().append('rect')
                .attr('class', 'ptm')
                .attr('id', function(d) { return d.value.residue + d.key; })
                .attr('x', function(d) { return axis(d.key - 1); })
                .attr('width', function(d) { return Math.max(ptm_viewer.min_barwidth, axis(1) - axis(0)); })
                .attr('y', function(d) { return ptm_viewer.yaxis( d.value.num_mods ); })
                .attr('height', function(d) { return ptm_viewer.yaxis( 0 ) - ptm_viewer.yaxis( d.value.num_mods ); })
                .attr('title', function(d) { return ptm_viewer.get_residue_tooltip_text(d); })
                .style('fill', function(d) { return residue_colors(d.value.residue); })
                .style('cursor', 'pointer')
                .on('mouseover', function() { d3.select(this).style('opacity', 0.8); })
                .on('mouseout', function() { d3.select(this).style('opacity', 1.0); })
                .on('click', function(d) { build_ptm_table(d.key, d.value, ptm_viewer.protein_data); });
};


PTMTrack.prototype.toggle_ptm = function(name, mode) {
    this.ptm_display_modes[name] = mode;
};

PTMTrack.prototype.toggle_exp = function(eid, mode) {
    this.experiment_display_modes[eid] = mode;
};

PTMTrack.prototype.get_residue_tooltip_text = function(data) {
    var site = data.key;
    var amino = data.value.residue;
    var mod_types = d3.keys(data.value.mods);
    var peptide = data.value.peptide;

    var text = "Residue: {0}{1}<br>15-mer: {2}<br>Modifications: {3}<br>".format(amino, site, peptide, mod_types.join(', '));

    return text;
}

PTMTrack.prototype.update_values = function(transition_duration) {
    ptm_viewer = this;
    t = this.g.transition().duration(transition_duration);

    t.selectAll('rect.ptm')
        .each(function(d){
            enabled = 0;
            mod_types = d3.keys(d.value.mods);
            for(var i in mod_types) {
                var k = mod_types[i];
                num_evidences = 0;
                for(var j in d.value.mods[k]){
                    exp_id = d.value.mods[k][j].experiment
                    if(ptm_viewer.experiment_display_modes[exp_id]) {
                        num_evidences+=1;
                    }
                }
                enabled += ptm_viewer.ptm_display_modes[k] && num_evidences > 0 ? 1 : 0;
            }
            d.value.num_mods = enabled;
        })
        .attr('y', function(d) { return ptm_viewer.yaxis(d.value.num_mods); })
        .attr('height', function(d) { return ptm_viewer.yaxis(0) - ptm_viewer.yaxis( d.value.num_mods ); });
}

PTMTrack.prototype.update_display = function(axis, viewer_width) {
    yaxis = this.yaxis;
    this.g.selectAll('rect.ptm')
        .attr('x', function(d) { return axis(d.key - 1); })
        .attr('width', function(d) { return Math.max(ptm_viewer.min_barwidth, axis(1) - axis(0));; })
        .attr('y', function(d) { return yaxis( d.value.num_mods ); })
        .attr('height', function(d) { return yaxis(0) - yaxis( d.value.num_mods ); });
};

function ResidueTrack(name, track_viewer, protein_data) {
    init_track(this, name, track_viewer, protein_data);

    // configurables
    this.height = 30;
    this.ticksize = 5;


    this.tick_levels = [5000,1000,500,100,50,10];
    ticks = [];
    for(var t in this.tick_levels)
        ticks.push([]);

    for(var i = 10; i < this.protein_data.seq.length; i+=10){
        for(var j in this.tick_levels){
            if(i % this.tick_levels[j] == 0){
                ticks[j].push(i);
                break;
            }
        }
    }
};

ResidueTrack.prototype.create = function(axis, viewer_width, show_residues) {
    this.show_residues = show_residues;

    this.g.append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', viewer_width)
            .attr('y1', 0)
            .attr('y2', 0);

    if(show_residues){
        this.g.selectAll('text.aminoacid')
                .data(this.protein_data.seq)
                    .enter().append('text')
                        .attr('class', 'aminoacid')
                        .attr('x', function(d,i){ return axis(i+0.5); })
                        .attr('y', "-0.1em")
                        .attr('text-anchor', 'middle')
                        .text(function(d) { return d; })
                        .style('pointer-events', 'none')
                        .style('font-size', axis(1) - axis(0));
    }

    for(var j in this.tick_levels){
        this.generate_ticks(viewer_width, axis, 5, ticks[j], this.tick_levels[j]);
    }

};

ResidueTrack.prototype.generate_ticks = function(viewer_width, axis, h, values, size){
    cls = 't{0}'.format(size);
    tick_opacity = 0;
    if(viewer_width / (axis(size) - axis(0)) < 20){
        tick_opacity = 1;
    }

    if(tick_opacity == 0)
        return;

    this.g.selectAll('line.'+cls)
        .data( values )
            .enter().append('line')
                .attr('class', cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); })
                .attr('y1', 0)
                .attr('y2', h)
                .style('opacity', tick_opacity);

    this.g.selectAll('text.'+cls)
        .data( values )
            .enter().append('text')
                .attr('class', cls)
                .attr('x', function(d) { return axis(d); })
                .attr('y', h)
                .attr('dy', '1em')
                .text(function(d) { return "" + d; } )
                .style('pointer-events', 'none')
                .style('opacity', tick_opacity);
};

ResidueTrack.prototype.update_display = function(axis, viewer_width) {
    for (var i in this.tick_levels){
        size = this.tick_levels[i];
        cls = 't{0}'.format(size);

        tick_opacity = 0;
        if(viewer_width / (axis(size) - axis(0)) < 20){
            tick_opacity = 1;
        }

        this.g.selectAll('line.'+cls)
                .attr('x1', function(d) { return axis(d); })
                .attr('x2', function(d) { return axis(d); })
                .style('opacity', tick_opacity);

        this.g.selectAll('text.'+cls)
                .attr('x', function(d) { return axis(d); })
                .style('opacity', tick_opacity);
    }

    if(this.show_residues){
        this.g.selectAll('text.aminoacid')
                    .attr('x', function(d,i){ return axis(i+0.5); })
                    .style('font-size', Math.min(16, axis(1) - axis(0)));
    }
};


function DomainTrack(name, track_viewer, protein_data) {
    init_track(this, name, track_viewer, protein_data);

    // configurables
    this.height = 40;
    this.domain_height = 20;

    this.g.append('text')
        .attr('x', 0)
        .attr('y', 0)
        .attr('dy', '-0.1em')
        .attr('text-anchor', 'left')
        .attr('class', 'track-label')
        .text(name);

};

DomainTrack.prototype.create = function(axis, viewer_width, domain_colors) {
    var pfam_url = this.protein_data.pfam_url;
    this.protein_domains = $.extend(true, [], this.protein_data.domains);

    this.g.append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', viewer_width)
            .attr('y1', 0)
            .attr('y2', 0);

    this.g.selectAll('rect.domain')
        .data(this.protein_domains)
            .enter().append('rect')
                .attr('class', 'domain')
                .attr('x', function(d) { return axis(d.start - 1); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start - 1); })
                .attr('y', 0)
                .attr('height', this.domain_height)
                .attr('title', function(d) { return d.label; })
                .style('fill', function(d) { return domain_colors( d.label ); } )
                .on('mouseover', function() { d3.select(this).style('opacity', 0.8); })
                .on('mouseout', function() { d3.select(this).style('opacity', 1.0); })
                .on('click', function(d) { window.open(pfam_url + d.label, '_blank'); });

    this.g.selectAll('text.domain')
        .data(this.protein_domains)
            .enter().append('text')
                .attr('class', 'domain')
                .attr('x', function(d) { return ( axis(d.start - 1) + axis(d.stop) ) / 2; })
                .attr('y', this.domain_height)
                .attr('dy', '1em')
                .attr('text-anchor', 'middle')
                .each(function(d) { d.show = ( axis(d.stop) - axis(d.start - 1) ) > d.label.length * 8; })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .text(function(d) { return d.label; });

};

DomainTrack.prototype.update_display = function(axis, viewer_width) {
    this.g.selectAll('rect.domain')
                .attr('x', function(d) { return axis(d.start - 1); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start - 1); });

    this.g.selectAll('text.domain')
                .each(function(d) { d.show = axis(d.stop) - axis(d.start - 1) > d.label.length * 8 })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .attr('x', function(d) { return (axis(d.start - 1) + axis(d.stop)) / 2; });
};


function EmptyTrack(name, track_viewer, protein_data) {
    init_track(this, name, track_viewer, protein_data);

    this.height = 50;
};

EmptyTrack.prototype.create = function(axis, viewer_width){

};

EmptyTrack.prototype.update_display = function(axis, viewer_width) {

};


function RegionTrack(name, track_viewer, protein_data) {
    init_track(this, name, track_viewer, protein_data);

    // configurables
    this.height = 40;
    this.region_height = 20;

};

RegionTrack.prototype.create = function(axis, viewer_width, region_colors) {
    var pfam_url = this.protein_data.pfam_url;
    this.protein_regions = $.extend(true, [], this.protein_data.regions);

    this.g.append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', viewer_width)
            .attr('y1', this.height)
            .attr('y2', this.height);

    this.g.selectAll('rect.region')
        .data(this.protein_regions)
            .enter().append('rect')
                .attr('class', 'region')
                .attr('x', function(d) { return axis(d.start-1); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start-1); })
                .attr('y', this.height - this.region_height)
                .attr('height', this.region_height)
                .attr('title', function(d) { return d.label; })
                .style('fill', function(d) { return region_colors( d.label ); } )
                .on('mouseover', function() { d3.select(this).style('opacity', 0.8); })
                .on('mouseout', function() { d3.select(this).style('opacity', 1.0); });

    this.g.selectAll('text.region')
        .data(this.protein_regions)
            .enter().append('text')
                .attr('class', 'region')
                .attr('x', function(d) { return ( axis(d.start-1) + axis(d.stop) ) / 2; })
                .attr('y', this.height - this.region_height)
                .attr('dy', "-0.25em")
                .attr('text-anchor', 'middle')
                .each(function(d) { d.show = ( axis(d.stop) - axis(d.start-1) ) > d.label.length * 8; })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .text(function(d) { return d.label; });
};

RegionTrack.prototype.update_display = function(axis, viewer_width) {
    this.g.selectAll('rect.region')
                .attr('x', function(d) { return axis(d.start-1); })
                .attr('width', function(d) { return axis(d.stop) - axis(d.start-1); });

    this.g.selectAll('text.region')
                .each(function(d) { d.show = axis(d.stop) - axis(d.start-1) > d.label.length * 8 })
                .style('opacity', function(d) { return d.show ? 1 : 0; })
                .attr('x', function(d) { return (axis(d.start-1) + axis(d.stop)) / 2; });
};

function MutationTrack(name, track_viewer, protein_data){
    init_track(this, name, track_viewer, protein_data);

    // configurables
    this.height = 45;
    this.min_mutation_size = 5;

    var max_mutations = 0;
    for(var k in protein_data.mutations) {
        var num_mutations = protein_data.mutations[k].length;
        protein_data.mutations[k].num_mutations = num_mutations;
        if(num_mutations > max_mutations){
            max_mutations = num_mutations;
        }
    }

    this.raxis = d3.scale.linear().domain([0, max_mutations]).range([5, 10]);

    this.g.append('text')
        .attr('x', 0)
        .attr('y', 0)
        .attr('dy', '0.5em')
        .attr('text-anchor', 'left')
        .attr('class', 'track-label')
        .text(name);

}

function get_mutation_tooltip(d){
    var tooltip_text = "Mutated Residue {0}{1}:<br>".format(d.value[0].original, d.key);

    var i = 0
    while(i < d.value.length){
        mutation = d.value[i];
        tooltip_text += "-> {0} | {1}<br>".format(mutation.mutant, mutation.annotation);

        i += 1;
    }

    return tooltip_text;
}

MutationTrack.prototype.create = function(axis, viewer_width, show_residues) {
    this.protein_regions = $.extend(true, [], this.protein_data.regions);
    var residue_max_width = (axis(1) - axis(0)) * 0.75;
    if(residue_max_width < this.min_mutation_size){
        residue_max_width = this.min_mutation_size;
    }
    this.raxis.range([residue_max_width/2, residue_max_width]);
    var raxis = this.raxis;

    this.show_residues = show_residues;

    this.g.append('line')
            .attr('class', "strand")
            .attr('x1', 0)
            .attr('x2', viewer_width)
            .attr('y1', this.height/2)
            .attr('y2', this.height/2);

    this.g.selectAll('circle.mutation')
        .data(d3.entries( this.protein_data.mutations ))
            .enter().append('circle')
                .attr('class', 'mutation')
                .attr('cx', function(d) { return axis(parseInt(d.key) - 0.5); })
                .attr('cy', this.height / 2 )
                .attr('r', function(d) { return raxis(d.value.length); })
                .attr('title', function(d) { return get_mutation_tooltip(d); })
                .style('fill', 'red' )
                .on('mouseover', function() { d3.select(this).style('fill', 'black'); })
                .on('mouseout', function() { d3.select(this).style('fill', 'red'); });

    if(this.show_residues){
        var original_residues = {};
        for(var k in this.protein_data.mutations){
            index = parseInt(k);
            original_residues[index] = this.protein_data.seq[index-1];
        }

        this.g.selectAll('text.mutation')
            .data(d3.entries(original_residues))
                .enter().append('text')
                    .attr('class', 'mutation')
                    .attr('x', function(d) { return axis(d.key - 0.5); })
                    .attr('y', this.height/2)
                    .attr('dy', '0.35em')
                    .attr('text-anchor', 'middle')
                    .style('pointer-events', 'none')
                    .style('fill', 'white')
                    .style('font-size', Math.min(16, axis(1) - axis(0)))
                    .text(function(d) { return d.value; });
    }
};

MutationTrack.prototype.update_display = function(axis, viewer_width) {
    var residue_max_width = (axis(1) - axis(0)) * 0.75;
    if(residue_max_width < this.min_mutation_size){
        residue_max_width = this.min_mutation_size;
    }
    this.raxis.range([residue_max_width/2, residue_max_width]);
    var raxis = this.raxis;

    this.g.selectAll('circle.mutation')
                .attr('cx', function(d) { return axis(parseInt(d.key) - 0.5); })
                .attr('r', function(d) { return raxis(d.value.length); });

    this.g.selectAll('text.mutation')
                .attr('x', function(d) { return axis(d.key - 0.5); })
                .style('font-size', Math.min(16, axis(1) - axis(0)));
};
