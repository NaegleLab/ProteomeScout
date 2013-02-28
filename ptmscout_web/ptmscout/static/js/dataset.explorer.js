function EnrichmentTable(parent_element, table_data){
	var widget = this;
	this.numTests = table_data.length;
	
	
	var div = $('<div />').appendTo(parent_element);
	$('<span>Mark as significant with p-value less than: </span>').appendTo(div);
	
	var last_val = "0.01";
	this.pvalueCutoff = $('<input />', {'value':last_val, 'type':"text",'length':"7"}).appendTo(div);
	this.pvalueCutoff.on('change', function(){
		var val = $(this).val();
		
		if(isNaN(parseFloat(val)) || !isFinite(val)){
			$(this).val(last_val);
		} else {
			last_val = val;
		}
	});
	this.setOptions = $('<button>Go</button>').appendTo(div);
	this.setOptions.on('click',function(){
		var cutoff_val = parseFloat(widget.pvalueCutoff.val());
		var correction = 'none'
		if(widget.bonferroniOption.is(":checked")){
			correction = 'bonferroni'
		}else if(widget.fdrOption.is(":checked")){
			correction = 'fdr'
		}
		
		
		widget.filter(cutoff_val, correction);
	});
	
	div = $('<div />').appendTo(parent_element);
	this.noCorrectionOption = $('<input />', {'type':"radio", 'name':"test_correction", 'value':"none", 'checked':"checked"}).appendTo(div);
	$('<span />', {'text':"No test correction"}).appendTo(div);
	
	div = $('<div />').appendTo(parent_element);
	this.bonferroniOption = $('<input />', {'type':"radio", 'name':"test_correction", 'value':"bonferoni"}).appendTo(div);
	$('<span />', {'text':"Bonferroni correction"}).appendTo(div);
	
	div = $('<div />').appendTo(parent_element);
	this.fdrOption = $('<input />', {'type':"radio", 'name':"test_correction", 'value':"fdr"}).appendTo(div);
	$('<span />', {'text':"False Discovery Rate correction"}).appendTo(div);
	
	this.tableElement = $('<table />', {'class':"enrichment-table"}).appendTo(parent_element);
	$('<thead><tr><th>Category</th><th>Term</th><th>P-Value</th></tr></thead>').appendTo(this.tableElement);
	this.tbody = $('<tbody />').appendTo(this.tableElement);
	
	for(var i in table_data){
		var value = table_data[i][2];
		var formatted_value = parseFloat(value.toPrecision(3)).toExponential()
		
		var row = $("<tr />").appendTo(this.tbody);
		$("<td />",{'text': table_data[i][0]}).appendTo(row);
		$("<td />",{'text': table_data[i][1]}).appendTo(row);
		$("<td />",{'value':value, 'text': formatted_value, 'class':"numeric"}).appendTo(row);
	}
	
	this.filter(parseFloat(last_val), 'none');
};

EnrichmentTable.prototype.filter = function(cutoff, correction) {
	var widget = this;
	
	this.tbody.find('tr').each(function() {
		var pval_td = $(this).find(".numeric");
		var val = pval_td.attr('value');
		var cval = val;
		
		
		if(correction == 'bonferroni')
			cval = val * widget.numTests;
		
		pval_td.text(parseFloat(parseFloat(cval).toPrecision(3)).toExponential());
		
		if(cval < cutoff)
			$(this).addClass('significant');
		else
			$(this).removeClass('significant');
		
		if(val > cutoff)
			$(this).addClass('notsignificant');
		else
			$(this).removeClass('notsignificant');
	});
};

function format_op_arg(arg){
	var oper_values = {'in':"IN", 'nin':"NOT IN", 'eq':"=", 'neq':"\u2260", 'gt':">", 'geq':"\u2265", 'lt':"<", 'leq':"\u2264", 'add':"+", 'sub':"-", 'mul':"x", 'div':"/"};

	if(arg in oper_values){
		return $('<span class="condition-op">{0}</span>'.format(oper_values[arg]));
	}
	return $('<span class="condition-arg">{0}</span>'.format(arg));
}

function format_query(query){
	if(query == 'experiment')
		return $('<div class="query-clause"><span class="condition-boolean">-</span><span class="condition-title">Entire Experiment</span></div>')
	
	var complete_query = $('<div></div>');
	
	for(var i in query){
		var bop = query[i][0];
		if(bop == 'nop'){
			bop = '-';
		}
		
		var condition = query[i][1];
		var title = condition[0];
		
		var qclause = $('<div class="query-clause"/>').appendTo(complete_query);
		
		$('<span class="condition-boolean">{0}</span>'.format(bop.toUpperCase())).appendTo(qclause);
		$('<span class="condition-title">{0}:</span>'.format(title.capitalize())).appendTo(qclause);
		
		for(var i = 1; i < condition.length; i++){
			format_op_arg( condition[i] ).appendTo( qclause );
		}
	}
	
	return complete_query;
}

function SubsetTab(tabId, name, selectionEngine) {
	var tabTemplate = "<li><a href='#{0}'>{1}</a> <span class='ui-icon ui-icon-close' role='presentation'>Remove Tab</span></li>";
	this.tabId = "tabs-" + tabId;
	
	var subsetTab = this;
	this.selectionEngine = selectionEngine;
	
	this.tabElement = $( tabTemplate.format(this.tabId, name) );
	this.contentContainer = $( '<div class="subset-content" id="{0}" />'.format(this.tabId) );
	
	this.filterset = $('<div class="filter-set" />').appendTo(this.contentContainer);
	this.fquery = $('<fieldset class="query"><legend>Foreground Query</legend></fieldset>').appendTo(this.filterset);
	this.bquery = $('<fieldset class="query"><legend>Background Query</legend></fieldset>').appendTo(this.filterset);
	var tmpdiv = $('<div class="formsubmission"></div>').appendTo(this.filterset);
	$('<div />', {'class':"clear"}).appendTo(this.filterset);
	
	$('<span>Subset Name:</span>').appendTo(tmpdiv);
	this.subsetName = $('<input type="text" value="{0}" />'.format(name)).appendTo(tmpdiv);
	this.saveSubset = $('<button></button>', {'class':'longtask', 'text':"Save Subset"}).appendTo(tmpdiv);
	this.saveSubset.on('click', function(){
		var name = subsetTab.subsetName.val();
		subsetTab.selectionEngine.saveSubset(subsetTab, name, subsetTab.response.foreground.query, subsetTab.response.background.query);
	});
	
	this.summary = $('<fieldset class="result-summary"><legend>Results</legend></fieldset>').appendTo(this.contentContainer);
	this.enrichment = $('<fieldset class="enrichment-tables"><legend>Feature Enrichment</legend></fieldset>').appendTo(this.contentContainer);
	
	this.frequency = $('<fieldset class="frequency-charts"><legend>Residue Frequency</legend><table><thead><tr><th>Foreground:</th><th>Background:</th></tr></thead><tbody><tr><td class="foreground-seqlogo"></td><td class="background-seqlogo"></td></tr></tbody></table></div>').appendTo(this.contentContainer);
	this.motif = $('<fieldset class="motifs"><legend>Motif Enrichment</legend></div>').appendTo(this.contentContainer);
	this.motifTable = $('<table><thead><tr><th>Motif</th><th>In Foreground</th><th>In Background</th><th>P-value</th></tr></thead></table>').appendTo(this.motif);
	this.motifEntries = $('<tbody></tbody>').appendTo(this.motifTable);

	this.proteins = $('<fieldset><legend>Proteins</legend></fieldset>').appendTo(this.contentContainer);
	this.proteinTable = $('<table class="protein_list"><thead><tr><th>MS id</th><th>Accession</th><th>Gene</th><th>Peptide</th><th>Site</th><th>Modifications</th></tr></thead></table>').appendTo(this.proteins);
	this.proteinEntries = $('<tbody></tbody>').appendTo(this.proteinTable);
	
	this.dataContainer = $('<fieldset><legend>Experiment Data</legend></fieldset>').appendTo(this.contentContainer);
	this.data = $('<div class="experiment_data"></div>').appendTo(this.dataContainer);
};

SubsetTab.prototype.changeName = function(name) {
	this.tabElement.find('a').text(name);
	this.subsetName.val(name);
};

SubsetTab.prototype.formatPeptides = function(peptides) {
	for(var i in peptides){
		var pepEntry = peptides[i];
		var rowClass = 'odd';
		if(i % 2 == 0)
			rowClass = 'even';
		
		var tableRow = $('<tr />', {"class":rowClass}).appendTo(this.proteinEntries);
		$('<td>{0}</td>'.format(pepEntry[0])).appendTo(tableRow);
		$('<td>{0}</td>'.format(pepEntry[1])).appendTo(tableRow);
		$('<td>{0}</td>'.format(pepEntry[2])).appendTo(tableRow);
		$('<td class="peptide">{0}</td>'.format(pepEntry[3])).appendTo(tableRow);
		$('<td class="modsite">{0}</td>'.format(pepEntry[4])).appendTo(tableRow);
		$('<td class="modsite">{0}</td>'.format(pepEntry[5])).appendTo(tableRow);
	}
};

SubsetTab.prototype.formatMotif = function(motif) {
	for(var i in motif.results){
		var rowClass = 'odd';
		if(i % 2 == 0)
			rowClass = 'even';
		var tableRow = $('<tr />', {"class":rowClass}).appendTo(this.motifEntries);
		var pval = parseFloat(motif.results[i][3].toPrecision(3)).toExponential();
		$('<td />',{ 'class':"peptide", 'text':motif.results[i][0] }).appendTo(tableRow);
		$('<td />',{ 'class':"numeric", 'text':'{0} / {1}'.format(motif.results[i][1][0], motif.results[i][1][1])}).appendTo(tableRow);
		$('<td />',{ 'class':"numeric", 'text':'{0} / {1}'.format(motif.results[i][2][0], motif.results[i][2][1])}).appendTo(tableRow);
		$('<td />',{ 'class':"numeric", 'text': pval}).appendTo(tableRow);
	}
};

SubsetTab.prototype.formatData = function(measurement_data){
	for(var i in measurement_data){
		var run_div = $('<div class="run" />').appendTo(this.data);
		var run = measurement_data[i];
		
		$('<div class="name">{0}</div>'.format(run.name)).appendTo(run_div);
		$('<div class="units">{0}</div>'.format(run.units)).appendTo(run_div);
		$('<div class="peptides">{0}</div>'.format(run.label)).appendTo(run_div);
		
		for(var j in run.series){
			var datapt = run.series[j];
			$('<div class="datapoint">{0},{1},{2}</div>'.format(datapt[0],datapt[1],datapt[2])).appendTo(run_div);
		}
	}
	
	$('<div class="chart" />').appendTo(this.data);
	d3.select("#{0}".format(this.tabId)).select(".experiment_data").each(function(){
		processRuns(this);
	});
}

SubsetTab.prototype.init = function(query_response) {
	format_query(query_response.foreground.query).appendTo(this.fquery);
	format_query(query_response.background.query).appendTo(this.bquery);
	
	this.response = query_response;
	
	$("<div>{0} out of {1} proteins selected</div>".format(query_response.foreground.proteins, query_response.background.proteins)).appendTo(this.summary)
	$("<div>{0} out of {1} peptides selected</div>".format(query_response.foreground.peptides, query_response.background.peptides)).appendTo(this.summary)
	$("<div>{0} out of {1} sites selected</div>".format(query_response.foreground.sites, query_response.background.sites)).appendTo(this.summary)

	this.enrichmentTable = new EnrichmentTable(this.enrichment, query_response.enrichment)
	
	createSeqlogo(d3.select("#{0}".format(this.tabId)).select(".foreground-seqlogo".format(this.tabId)), query_response.foreground.seqlogo, 500, 400);
	createSeqlogo(d3.select("#{0}".format(this.tabId)).select(".background-seqlogo".format(this.tabId)), query_response.background.seqlogo, 500, 400);
	
	this.formatMotif(query_response.motif);
	this.formatPeptides(query_response.peptides);
	this.formatData(query_response.measurements);
};

function QuantitativeSelector(parent_element, field_values) {
	this.element = $("<span>:</span>");
	this.element.appendTo(parent_element);

	this.field_values = field_values;
	this.compare_values = {'eq':"=", 'neq':"\u2260", 'gt':">", 'geq':"\u2265", 'lt':"<", 'leq':"\u2264"};
	this.oper_values = {'add':"+", 'sub':"-", 'mul':"x", 'div':"/"}
	
	this.selectors = [];
	this.oper_state = 'LHS';
	
	this.addVariable();
};

QuantitativeSelector.prototype.build_query = function(){
	if(!this.checkValid()){
		throw new Error("Invalid expression");
	}
	
	var expression = [];
	expression.push("quantitative");
	
	for(var i in this.selectors){
		if(this.selectors[i].hasClass['operator']){
			var v = this.selectors[i].val();
			expression.push(v);
		}
		if(this.selectors[i].hasClass['variable']){
			var v = this.selectors[i].find('select').val();
			var n = this.selectors[i].find('input').val();
			
			if(v == ''){
				throw new Error("Variable cannot be empty");
			}
			
			if(v == 'numeric'){
				if(n == ''){
					throw new Error("Numeric value cannot be empty");
				}
				expression.push(n);
			}else
				expression.push(v);
		}
	}
	
	var state = this.computeOperState(this.selectors.length-1);
	if(state == 'L'){
		throw new Error("Incomplete quantiative comparison");
	}
	
	if(expression[expression.length-1] == ''){
		expression.pop();
	}
	
	if(expression.length == 0){
		throw new Error("Empty Expression");
	}

	return expression;
};

QuantitativeSelector.prototype.computeOperState = function(oper_index) {
	var cstate = ''; 
	for(var i = 1; i <= oper_index; i+=2){
		var oper_val = this.selectors[i].val();

		if(cstate == 'R'){
			cstate = 'D';
		}
		if(cstate == 'M'){
			cstate = 'R';
		}
		if(cstate == 'L'){
			cstate = 'M';
		}
		if(cstate == ''){
			cstate = 'L';
		}
		
		if(oper_val in this.compare_values){
			cstate = 'M';
		}
	}
	return cstate;
};

QuantitativeSelector.prototype.checkValid = function() {
	var last_val = "start";
	for(var i in this.selectors){
		if(this.selectors[i].hasClass('variable')){
			var var_val = this.selectors[i].find('select').val();
			if(last_val == ''){
				return false;
			}
			if(var_val == ''){
				return false;
			}
		}
		if(this.selectors[i].hasClass('operator')){
			var oper_state = this.computeOperState(i);
			var oper_val = this.selectors[i].val();
			
			if(last_val == ''){
				return false;
			}
			if(oper_state == 'D'){
				return false;
			}
			if(oper_state == 'M' && !(oper_val in this.compare_values) && oper_val != ''){
				return false;
			}
			if(oper_state == 'R' && !(oper_val in this.oper_values) && oper_val != ''){
				return false;
			}
			
			last_val = oper_val;
		}
	}
	return true;
};

QuantitativeSelector.prototype.removeAfter = function(i) {
	while(this.selectors.length > i+1){
		var elem = this.selectors.pop();
		elem.remove();
	}
};

QuantitativeSelector.prototype.addOperator = function() {
	var selector = this;
	
	var i = this.selectors.length;

	var oper_select = $("<select />", {class:"operator"});
	oper_select.on('change', function(){
		var val = $(this).val();
		
		if(!selector.checkValid()){
			selector.removeAfter(i);
		}
		if (i+1 == selector.selectors.length){
			selector.addVariable();
		}
	});
	this.selectors.push(oper_select);
	var oper_state = this.computeOperState(i+1);
	if(oper_state == 'D'){
		this.selectors.pop();
		return;
	}
	
	
	$("<option />").appendTo(oper_select);
	
	if(oper_state != 'M'){
		for(var v in this.oper_values){
			$("<option value=\"{0}\">{1}</option>".format( v, this.oper_values[v] )).appendTo(oper_select);
		}	
	}
	if(oper_state != 'R'){
		for(var v in this.compare_values){
			$("<option value=\"{0}\">{1}</option>".format( v, this.compare_values[v] )).appendTo(oper_select);
		}
	}	
	
	oper_select.appendTo(this.element);
};

QuantitativeSelector.prototype.addVariable = function() {
	var selector = this;
	
	var i = this.selectors.length;
	
	var container = $("<span />", {class:"variable"}); 
	var variable_select = $("<select />").appendTo(container);
	var numeric_input = $('<input type="text" style="display:none;"/>').appendTo(container);
	var last_val = "";
	
	numeric_input.on('change', function(){
		var val = $(this).val();
		
		if(isNaN(parseFloat(val)) || !isFinite(val)){
			$(this).val(last_val);
		} else {
			last_val = val;
		}
	});
	
	this.selectors.push(container);
	
	variable_select.on('change', function(){
		if($(this).val() == ''){
			selector.removeAfter(i);
		} 
		if($(this).val() == 'numeric'){
			numeric_input.show();
		} else {
			numeric_input.hide();
		}
		if(i+1 == selector.selectors.length && selector.computeOperState(i) != 'D'){
			selector.addOperator();
		}
		
	});
	
	$("<option />").appendTo(variable_select);
	$('<option value="numeric">Numeric</option>').appendTo(variable_select);
	for(var j in this.field_values){
		$("<option value=\"{0}\">{0}</option>".format( this.field_values[j] )).appendTo(variable_select);
	}
	
	container.appendTo(this.element);
};

function SubsetSelector(parent_element, saved_subsets) {
	this.element = $("<span />");
	
	this.element.appendTo(parent_element);
	
	$("<span>MS ID</span>").appendTo(this.element);
	
	this.operation = $("<select />").appendTo(this.element);
	$("<option value=\"in\">in</option>").appendTo(this.operation);
	$("<option value=\"nin\">not in</option>").appendTo(this.operation);
	
	$("<span class=\"field_label\">Subset: </span>").appendTo(this.element);
	this.value = $("<select />").appendTo(this.element);
	
	$("<option />").appendTo(this.value);
	for(var i in saved_subsets){
		$("<option value=\"{0}\">{0}</option>".format(saved_subsets[i])).appendTo(this.value);
	}
};

SubsetSelector.prototype.build_query = function(){
	var op = this.operation.val();
	var subset = this.value.val()
	
	if(subset == '')
		throw new Error("Subset name must be specified");
	
	return ['subset', op, subset];
};

function ScansiteSelector(parent_element, scansite_keys, scansite_fields) {
	var selector = this;
	this.element = $("<span />");
	
	this.field_values = scansite_fields;
	this.element.appendTo(parent_element);
	
	$("<span>Type: </span>").appendTo(this.element);
	this.field = $("<select />").appendTo(this.element);
	
	$('<option />').appendTo(this.field);
	for(var i in scansite_keys){
		var k = scansite_keys[i];
		$('<option value="{0}">{0}</option>'.format(k)).appendTo(this.field);
	}
	this.field.on('change', function(){
		var field_val = selector.field.val();
		selector.value.empty();
		
		if(field_val in selector.field_values){
			$("<option />").appendTo(selector.value);
			for(var i in selector.field_values[field_val]){
				$("<option value=\"{0}\">{0}</option>".format( selector.field_values[field_val][i] )).appendTo(selector.value);
			}
		}
	});
	
	
	this.operation = $("<select />").appendTo(this.element);
	$("<option value=\"eq\">=</option>").appendTo(this.operation);
	$("<option value=\"neq\">\u2260</option>").appendTo(this.operation);
	
	$("<span class=\"field_label\">Term: </span>").appendTo(this.element);
	this.value = $("<select />").appendTo(this.element);
//	this.value.combobox();
	
	$("<span class=\"field_label\">Stringency: </span>").appendTo(this.element);
	this.stringency = $("<select />").appendTo(this.element);
	
	$("<option value=\"5\">Low</option>").appendTo(this.stringency);
	$("<option value=\"1\">Medium</option>").appendTo(this.stringency);
	$("<option value=\"0.2\">High</option>").appendTo(this.stringency);
};

ScansiteSelector.prototype.build_query = function() {
	var f = this.field.val();
	var op = this.operation.val();
	var term = this.value.val();
	var s = this.stringency.val();
	
	if(f == '')
		throw new Error("You must select a type");
	if(term == '')
		throw new Error("You must select a term");
	
	return ['scansite', f, op, term, s];
};

function ProteinSelector(parent_element, field_values) {
	var selector = this;
	
	this.element = $("<span />");
	this.element.appendTo(parent_element);
	$("<span class=\"field_label\">Gene/Protein: </span>").appendTo(this.element);
	this.value = $("<input type=\"text\" width=\"200\" />").appendTo(this.element);
	
	this.value.autocomplete({
	      source: field_values
	    });
	
//	$('<option />').appendTo(this.value);
//	for(var i in field_values){
//		$("<option value=\"{0}\">{0}</option>".format(field_values[i])).appendTo(this.value);
//	}
//	this.value.combobox();
};

ProteinSelector.prototype.build_query = function() {
	var p = this.value.val();
	
	if(p == ''){
		throw new Error("Gene/Protein search field cannot be empty");
	}
	return ['protein', p];
};


function MetadataSelector(parent_element, type, field_name, value_name, field_values, values_by_field) {
	var selector = this;

	this.type = type;
	this.field_values = field_values;
	this.values_by_field = values_by_field;
	
	this.element = $("<span />");
	this.element.appendTo(parent_element);
	$("<span class=\"field_label\">{0}</span>".format(field_name)).appendTo(this.element);
	this.field_name = $("<select />").appendTo(this.element);
	
	$("<option />").appendTo(this.field_name);
	for(var i in field_values){
		$("<option value=\"{0}\">{0}</option>".format( field_values[i] )).appendTo(this.field_name);
	}
	
	this.operation = $("<select />").appendTo(this.element);
	$("<option value=\"eq\">=</option>").appendTo(this.operation);
	$("<option value=\"neq\">\u2260</option>").appendTo(this.operation);
	
	$("<span class=\"field_label\">{0}</span>".format(value_name)).appendTo(this.element);
	this.value = $("<select />").appendTo(this.element);
	
	this.field_name.on('change', function(){
		var field_val = selector.field_name.val();
		selector.value.empty();
		
		if(field_val in selector.values_by_field){
			$("<option />").appendTo(selector.value);
			for(var i in selector.values_by_field[field_val]){
				$("<option value=\"{0}\">{0}</option>".format( selector.values_by_field[field_val][i] )).appendTo(selector.value);
			}
		}
	});
//	this.value.combobox();
};

MetadataSelector.prototype.build_query = function() {
	var f = this.field_name.val();
	var op = this.operation.val();
	var v = this.value.val();
	
	if(f == '')
		throw new Error("You must specify a field");
	if(v == '')
		throw new Error("You must specify a value");
	
	return [this.type, f, op, v]
};

function SequenceSelector(parent_element) {
	this.element = $("<span />");
	this.element.appendTo(parent_element);
	
	$("<span class=\"field_label\">Peptide Sequence: </span>").appendTo(this.element);
	this.value = $("<input type=\"text\" length=\"15\" width=\"15\" />").appendTo(this.element);
};

SequenceSelector.prototype.build_query = function() {
	var v = this.value.val();
	if(v==''){
		throw new Error("You must specify a search string");
	}
	
	return ['sequence', v]
};

function SelectorCondition(parent_element, field_data, show_boolean_op) {
	var condition = this;
	this.field_data = field_data;
	
	this.element = $('<div></div>');
	var boolean_container = $('<span class="condition-boolean"></span>');
	boolean_container.appendTo(this.element);
	
	if(show_boolean_op){
		this.boolean_operator = $('<select><option value="and">and</option><option value="or">or</option></select>').appendTo(boolean_container);
	}else{
		this.boolean_operator = null;
		boolean_container.text("-");
	}
	$("<span class=\"condition-title\">Condition:</span>").appendTo(this.element);
	this.type_select = $("<select />").appendTo(this.element);
	
	$("<option></option>").appendTo(this.type_select);
	$("<option value=\"quantitative\">Quantitative</option>").appendTo(this.type_select);
	$("<option value=\"subset\">Subset</option>").appendTo(this.type_select);
	$("<option value=\"cluster\">Cluster</option>").appendTo(this.type_select);
	$("<option value=\"metadata\">Metadata</option>").appendTo(this.type_select);
	$("<option value=\"protein\">Protein</option>").appendTo(this.type_select);
	$("<option value=\"scansite\">Scansite</option>").appendTo(this.type_select);
	$("<option value=\"sequence\">Sequence</option>").appendTo(this.type_select);
	
	
	this.type_select.on('change', function(){
		condition.changed();
	});
	
	this.selector = null;
	
	this.element.appendTo(parent_element);
};

SelectorCondition.prototype.build_query = function() {
	var expression = [];
	if(this.boolean_operator != null)
		expression.push( this.boolean_operator.val() );
	else
		expression.push( 'nop' );
	
	var value = this.type_select.val();
	if(value == '')
		throw new Error("Empty condition clause");
	
	expression.push( this.selector.build_query() );
	
	return expression;
};

SelectorCondition.prototype.changed = function() {
	var value = this.type_select.val();
	if(this.selector != null){
		this.selector.element.remove();
	}
	if(value == 'quantitative'){
		this.selector = new QuantitativeSelector(this.element, this.field_data.quantitative_fields);
	}
	if(value == 'metadata'){
		this.selector = new MetadataSelector(this.element, 'metadata', ":", "", this.field_data.metadata_keys, this.field_data.metadata_fields);
	}
	if(value == 'subset'){
		this.selector = new SubsetSelector(this.element, this.field_data.subset_labels);
	}
	if(value == 'cluster'){
		this.selector = new MetadataSelector(this.element, 'cluster', ":", "Cluster ID:", this.field_data.clustering_sets, this.field_data.clustering_labels);
	}
	if(value == 'sequence'){
		this.selector = new SequenceSelector(this.element);
	}
	if(value == 'protein'){
		this.selector = new ProteinSelector(this.element, this.field_data.accessions);
	}
	if(value == 'scansite'){
		this.selector = new ScansiteSelector(this.element, this.field_data.scansite_keys, this.field_data.scansite_fields);
	}
};

function SubsetManager(element, selectionEngine) {
	this.tabCounter = 0;
	this.tabs = element.tabs();
	var tabs = this.tabs
	
	this.selectionEngine = selectionEngine;
	
	tabs.delegate( "span.ui-icon-close", "click", function() {
	      var panelId = $( this ).closest( "li" ).remove().attr( "aria-controls" );
	      $( "#" + panelId ).remove();
	      tabs.tabs( "refresh" );
	    });
	 
    tabs.bind( "keyup", function( event ) {
      if ( event.altKey && event.keyCode === $.ui.keyCode.BACKSPACE ) {
        var panelId = tabs.find( ".ui-tabs-active" ).remove().attr( "aria-controls" );
        $( "#" + panelId ).remove();
        tabs.tabs( "refresh" );
      }
    });
}

SubsetManager.prototype.addSubset = function(query_result) {
	var st = new SubsetTab(this.tabCounter, query_result.name, this.selectionEngine)
	
	var numtabs = this.tabs.find( ".ui-tabs-nav > li" ).length;
	this.tabs.find( ".ui-tabs-nav" ).append( st.tabElement );
    this.tabs.append( st.contentContainer );
    this.tabs.tabs( "refresh" );
    
    st.init(query_result);
    this.tabs.tabs('option', 'active', numtabs);
    
	this.tabCounter++;
};

function SubsetSelection(element, webservice_url) {
	var selector = this;
	
	this.element = element;
	
	this.webservice_url = webservice_url;
	this.query_num = 0;
	
	this.subsetPrefix = 'Subset: ';
	this.field_data = JSON.parse( Base64.decode( element.find("#field-data").text() ) );
	this.experiment_id = this.field_data.experiment_id;
	

	this.failure_dialog = $("<div title=\"Error!\"></div>")
								.dialog( {  autoOpen:false, modal: true, buttons: { "Ok": function(){ $(this).dialog( 'close' ); } } } );
	this.failure_text = $("<span />").appendTo(this.failure_dialog);
	
	this.subsetSelection = $("#saved-subset-select");
	this.backgroundSelection = $("#background-select");
	
	for(var i in this.field_data.subset_labels){
		var label = this.field_data.subset_labels[i];
		$('<option />', {'value': label, 'text': label}).appendTo(this.subsetSelection);
		var background_label  = '{0}{1}'.format(this.subsetPrefix, label);
		$('<option />', {'value': background_label, 'text': background_label}).appendTo(this.backgroundSelection);
	}
	
	this.subsetOpen = $("#open-saved-subset");
	this.subsetOpen.on('click', function(){
		selector.openExistingSubset(selector.subsetSelection.val());
	})
	
	
	this.conditions = [];
	this.condition_list = $("#filter-conditions");
	this.add_condition = $("#add-condition");
	this.add_condition.on('click', function(){
		selector.addConditionField();
	});
	this.remove_condition = $("#remove-condition");
	this.remove_condition.on('click', function(){
		selector.removeConditionField();
	});
	

	
	this.compute_subset = $("#compute-subset");
	this.compute_subset.on('click',function(){
		selector.submitQuery();
	});
	
	this.clear_form = $("#clear-form");
	this.clear_form.on('click', function(){
		selector.clearConditions();
	});
	
	this.subsetManager = new SubsetManager( element.find('#open-subsets'), this );
};

SubsetSelection.prototype.displayError = function(errorMessage) {
	this.failure_text.text(errorMessage);
	this.failure_dialog.dialog( "open" );
	done_waiting();
};

SubsetSelection.prototype.openExistingSubset = function(name) {
	var app = this;
	
	if(name == ''){
		this.displayError("You must choose a subset");
	}else{
		var query_url = this.webservice_url + '/subsets/fetch';
		var subset_fetch_query = {	'name': name,
									'experiment': this.experiment_id 
				 					}
		
		$.ajax({
			url: query_url,
			type: 'POST',
			contentType:'application/json',
			data: JSON.stringify(subset_fetch_query),
			dataType:'json',
			success: function(data) {
				done_waiting();
				
				if(data.status == 'error'){
					app.displayError(data.message);
				}else if(data.status == 'success'){
					app.subsetManager.addSubset(data);
				}
			} ,
			error: function(jqXHR, textStatus, errorThrown) {
				app.displayError("Error in server response, please try again.");
				console.log(textStatus);
				console.log(errorThrown);
			}
		});

	}
}

SubsetSelection.prototype.saveSubset = function(tabElement, name, foreground, background){
	var app = this;
	
	if(name == ''){
		this.displayError("Name cannot be an empty field");
	}else{

		var query_url = this.webservice_url + '/subsets/save';
		var subset_save_query = {'name': name,
								 'experiment': this.experiment_id,
								 'foreground': foreground,
								 'background': background}		
		
		$.ajax({
			url: query_url,
			type: 'POST',
			contentType:'application/json',
			data: JSON.stringify(subset_save_query),
			dataType:'json',
			success: function(data) {
				done_waiting();
				
				if(data.status == 'error'){
					app.displayError(data.message);
				}else if(data.status == 'success'){
					var label = data.name;
					tabElement.changeName(label);

					app.field_data.subset_labels.push(label);
					$('<option />', {'value': label, 'text': label}).appendTo(app.subsetSelection);
					
					var background_label  = 'Subset: {0}'.format(label);
					$('<option />', {'value': background_label, 'text': background_label}).appendTo(app.backgroundSelection);
					
				}
			} ,
			error: function(jqXHR, textStatus, errorThrown) {
				app.displayError("Error in server response, please try again.");
				console.log(textStatus);
				console.log(errorThrown);
			}
		});
		
	}
};

SubsetSelection.prototype.submitQuery = function() {
	var app = this;
	
	this.query_num += 1;
	if(this.conditions.length == 0){
		this.displayError("You must specify query conditions!");
		return;
	}
	
	var expression = [];
	for(var i in this.conditions){
		try{
			q = this.conditions[i].build_query()
			expression.push(q);
		}catch(err){
			this.displayError("Query error line {0}: {1}".format(parseInt(i)+1, err.message));
			return;
		}
	}
	
	var background = this.backgroundSelection.val();
	if(background.startsWith(this.subsetPrefix)){
		var bg_subset = background.substring(this.subsetPrefix.length);
		background = [ ['nop', ['subset', 'in', bg_subset]] ];
	}
	
	var subset_query = {
			'experiment': this.experiment_id,
			'type': 'create',
			'name': 'Subset {0}'.format(this.query_num),
			'background': background,
			'foreground': expression
	};
	
	var query_url = this.webservice_url + '/subsets/query';
	
	$.ajax({
		url: query_url,
		type: 'POST',
		contentType:'application/json',
		data: JSON.stringify(subset_query),
		dataType:'json',
		success: function(data) {
			app.subsetManager.addSubset(data);
			done_waiting();
		} ,
		error: function(jqXHR, textStatus, errorThrown) {
			app.displayError("Error in server response, please try again.");
			console.log(textStatus);
			console.log(errorThrown);
		}
	});
}

SubsetSelection.prototype.addConditionField = function() {
	var is_first = this.conditions.length == 0;
	var new_field = new SelectorCondition(this.condition_list, this.field_data, !is_first);
	this.conditions.push(new_field);
};

SubsetSelection.prototype.removeConditionField = function() {
	if(this.conditions.length > 0){
		var field = this.conditions.pop();
		field.element.remove();
	}
};

SubsetSelection.prototype.clearConditions = function() {
	for(var i in this.conditions){
		this.conditions[i].remove();
	}
	this.conditions = [];
};

$(function(){
	new SubsetSelection($("#subset-select"), $("#webservice_url").text());
});
