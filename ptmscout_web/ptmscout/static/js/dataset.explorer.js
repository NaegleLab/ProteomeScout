function QuantitativeSelector(parent_element, field_values) {
	this.element = $("<span>:</span>");
	this.element.appendTo(parent_element);

	this.field_values = field_values;
	this.compare_values = {'eq':"=", 'neq':"\u2260", 'gt':">", 'geq':"\u2265", 'lt':"<", 'leq':"\u2264"};
	this.oper_values = {'plus':"+", 'minus':"-", 'multiply':"x", 'divide':"/"}
	
	this.selectors = [];
	this.oper_state = 'LHS';
	
	this.addVariable();
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
	for(var i in this.selectors){
		if(this.selectors[i].hasClass('operator')){
			var oper_state = this.computeOperState(i);
			var oper_val = this.selectors[i].val();
			
			if(oper_state == 'D'){
				return false;
			}
			if(oper_state == 'M' && !(oper_val in this.compare_values)){
				return false;
			}
			if(oper_state == 'R' && !(oper_val in this.oper_values)){
				return false;
			}
		}
	}
	return true;
}

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
	
	var variable_select = $("<select />", {class:"variable"});
	this.selectors.push(variable_select);
	
	variable_select.on('change', function(){
		if($(this).val() == ''){
			selector.removeAfter(i);
		} 
		if(i+1 == selector.selectors.length && selector.computeOperState(i) != 'D'){
			selector.addOperator();
		}
		
	});
	
	$("<option />").appendTo(variable_select);
	for(var j in this.field_values){
		$("<option value=\"{0}\">{0}</option>".format( this.field_values[j] )).appendTo(variable_select);
	}
	
	variable_select.appendTo(this.element);
};

function SubsetSelector(parent_element) {
	this.element = $("<span />");
	
	this.element.appendTo(parent_element);
	
	this.operation = $("<select />").appendTo(this.element);
	$("<option value=\"in\">in</option>").appendTo(this.operation);
	$("<option value=\"not in\">not in</option>").appendTo(this.operation);
	
	$("<span class=\"field_label\">Subset: </span>").appendTo(this.element);
	this.value = $("<select />").appendTo(this.element);
};

function MetadataSelector(parent_element, field_name, value_name, field_values, values_by_field) {
	var selector = this;

	this.field_values = field_values;
	this.values_by_field = values_by_field;
	
	this.element = $("<span />");
	this.element.appendTo(parent_element);
	$("<span class=\"field_label\">{0}: </span>".format(field_name)).appendTo(this.element);
	this.field_name = $("<select />").appendTo(this.element);
	
	for(var i in field_values){
		$("<option value=\"{0}\">{0}</option>".format( field_values[i] ));
	}
	

	
	this.operation = $("<select />").appendTo(this.element);
	$("<option value=\"in\">in</option>").appendTo(this.operation);
	$("<option value=\"not in\">not in</option>").appendTo(this.operation);
	
	$("<span class=\"field_label\">{0}: </span>".format(value_name)).appendTo(this.element);
	this.value = $("<select />").appendTo(this.element);
	
	this.field_name.on('change', function(){
		selector.value.empty();
		for(var i in values_by_field){
			$("<option value=\"{0}\">{0}</option>".format( selector.values_by_field[i] ));
		}
	});
};

function SequenceSelector(parent_element) {
	this.element = $("<span />");
	this.element.appendTo(parent_element);
	
	$("<span class=\"field_label\">Peptide Sequence: </span>").appendTo(this.element);
	this.value = $("<input type=\"text\" length=\"15\" width=\"15\" />").appendTo(this.element);
};

function SelectorCondition(parent_element, field_data) {
	var condition = this;
	
	this.field_data = field_data;
	
	this.element = $('<div><span class=\"condition-title\">Condition:</span></div>');
	this.type_select = $("<select />").appendTo(this.element);
	
	$("<option></option>").appendTo(this.type_select);
	$("<option value=\"quantitative\">Quantitative</option>").appendTo(this.type_select);
	$("<option value=\"subset\">Subset</option>").appendTo(this.type_select);
	$("<option value=\"cluster\">Cluster</option>").appendTo(this.type_select);
	$("<option value=\"metadata\">Metadata</option>").appendTo(this.type_select);
	$("<option value=\"sequence\">Sequence</option>").appendTo(this.type_select);
	
	this.type_select.on('change', function(){
		condition.changed();
	});
	
	this.selector = null;
	
	this.element.appendTo(parent_element);
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
		this.selector = new MetadataSelector(this.element, "Metadata Field", "", this.field_data.metadata_keys, this.field_data.metadata_fields);
	}
	if(value == 'subset'){
		this.selector = new SubsetSelector(this.element);
	}
	if(value == 'cluster'){
		this.selector = new MetadataSelector(this.element, "Cluster Set", "Cluster ID", this.field_data.clustering_sets, this.field_data.clustering_labels);
	}
	if(value == 'sequence'){
		this.selector = new SequenceSelector(this.element);
	}
};

SelectorCondition.prototype.remove = function() {
	this.element.remove();
};

function SubsetSelection(element) {
	var selector = this;
	this.element = element;
	
	this.field_data = JSON.parse( Base64.decode( element.find("#field-data").text() ) );
	
	this.conditions = [];
	this.condition_list = $("#filter-conditions");
	this.add_condition = $("#add-condition");
	this.add_condition.on('click', function(){
		selector.add_condition_field();
	});
	
	this.compute_subset = $("#compute-subset");
	this.clear_form = $("#clear-form");
	this.clear_form.on('click', function(){
		selector.clear_conditions();
	});
	
	this.open_subsets = element.find('#open-subsets').tabs();
};

SubsetSelection.prototype.add_condition_field = function() {
	var new_field = new SelectorCondition(this.condition_list, this.field_data);
	this.conditions.push(new_field);
};

SubsetSelection.prototype.clear_conditions = function() {
	for(var i in this.conditions){
		this.conditions[i].remove();
	}
	this.conditions = [];
};

$(function(){
	new SubsetSelection($("#subset-select"));
});
