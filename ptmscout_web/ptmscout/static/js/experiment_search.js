function ExperimentEntry(table_element, edata) {
    this.marked = false;
    this.data = edata;
    this.element = $("<tr></tr>", {
                        'class': "entry",
                        'id': "e{0}".format(this.data.id)
                    });

    table_element.append(this.element);
    $(this.element).on('click', null, this, function(e) {
                                e.data.mark_selected();
                            });

    $("<td><input type=\"hidden\" value=\"e{0}\" name=\"experiment\" /></td>".format(edata.id)).appendTo(this.element);
    $("<td>{0}</td>".format(edata.id)).appendTo(this.element);
    $("<td>{0}</td>".format(edata.name)).appendTo(this.element);
    $("<td>{0}</td>".format(edata.residues)).appendTo(this.element);
};

ExperimentEntry.prototype.remove = function() {
    this.element.remove();
};

ExperimentEntry.prototype.mark_selected = function() {
    if(!this.marked){
        this.marked = true;
        this.element.addClass('marked');
    }else{
        this.marked = false;
        this.element.removeClass('marked');
    }
};

function ExperimentSearch(finder_element, list_element, service_url){
    var search_form = this;
    this.finder = finder_element;
    this.list = list_element;
    this.results = this.finder.find("#result_list");

    this.added_list = this.list.find(".experiment_list > tbody");
    this.result_list = this.finder.find("#result_list > table > tbody");
    this.text_search = this.finder.find("#text_search");
    this.submit_button = this.finder.find("#submit_query");

    this.add_conditions = this.finder.find("#add_conditions");
    this.clear_conditions = this.finder.find("#clear_conditions");

    this.add_button = this.finder.find("#add_selected");

    this.add_button.on( 'click', function() {
        search_form.add_selection();
    });

    this.submit_button.on( 'click', function() {
        search_form.submit_search(service_url);
    });

    this.search_results = [];
    this.added_results = [];
};

ExperimentSearch.prototype.submit_search = function(service_url) {
    var search_form = this;
    var search_value = this.text_search.val();

    var url_args = { 'search_term': search_value };

    console.log( url_args );


    search_form.result_list.find(".entry").remove();
    this.search_results = [];

    $.ajax({
          dataType: "json",
          url: "{0}/{1}".format(service_url, "experiments"),
          data: url_args,
          success: function(data){
              for(var i in data.experiments){
                  search_form.search_results.push(new ExperimentEntry(search_form.result_list, data.experiments[i]));
              }
              search_form.results.show()
              $("#waiting-dialog").dialog("close");
          }
    });
};

ExperimentSearch.prototype.add_selection = function() {
    var results = this.search_results;

    for(var i in results){
        if(results[i].marked)
            this.added_results.push(new ExperimentEntry(this.added_list, results[i].data));
    }

    if(this.added_results.length > 0){
        this.list.show();
    }
};

ExperimentSearch.prototype.remove_selection = function() {
    var nlist = [];
    for(var i in this.added_results){
        if(this.added_results[i].marked){
            this.added_results[i].remove();
        }
        else
            nlist.push(this.added_results[i]);
    }
    this.added_results = nlist;
    if(this.added_results.length == 0){
        this.list.hide();
    }
};

ExperimentSearch.prototype.show = function() {
    var esearch = this;
    this.finder.dialog({
                width: 800,
                height: 500
                });
};

$(function() {
    $("#chosen_list").hide();
    $("#result_list").hide();
    var web_url = $("#webservice_url").text();

    var esearch = new ExperimentSearch($("#experiment_search"), $("#chosen_list"), web_url);

    $("#choose_experiments").click(function(){
        esearch.show();
    });
    $("#remove_selected").click(function(){
        esearch.remove_selection();
    });
});
