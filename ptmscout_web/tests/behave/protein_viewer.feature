Feature: Protein Viewer
	 In order to easily explore proteins
	 Users should be able to manipulate the view of a protein


	Scenario: Zoom in on a protein
		  Given a user clicks on the magnify in function
		  When the user clicks on a protein area or one of the protein track areas
		  Then the lateral portion of the protein and its associated tracks should expand to be centered on the zone of the click and doubling in zoom size with each click

	Scenario: Zoom limit
		  Given a user clicks on a magnify function and the protein is either all the way zoomed out or all the way zoomed in
		  When the user clicks on the protein area
		  Then do nothing

	Scenario: Zoom out on a protein 
		  Given a user clicks on the magnify out function
		  When the user clicks on a protein area or one of the protein track areas
		  Then the lateral portion of the protein and protein tracks should contract to be centered on the zone of the click and halving in zoom size with each click

	Scenario: Mouse-over modifications
		  Given a user is viewing a protein with modifications
		  When the user mouses over a particular site
		  Then a descriptive box should appear with the known modifications on that site, the 15-mer sequence centered on that site, and with a link to the experimental evidence of that site's modifications

   	Scenario: Mouse-over click
		  Given a user has moused over a site
		  When the user clicks on a modification of that site
		  Then the user should be taken to the portion of the experimental table detailing the evidences of that modification

	Scenario: Filter by modification
		  Given a user is looking at a protein with many types of modifications
		  When the user de-selects a modification type 
		  Then the indication of that modification appearing on sites within the protein should disappear
		
	Scenario: Turn on or off tracks
		  Given a user is looking at a protein with multiple track types
		  When the user toggles the track box
		  Then the track should appear or disappear accordingly

	Scenario: Export peptide sequence      #where to put this feature
		  Given a user requests the export of a peptide sequence
		  When some modifications have been filtered
		  Then export the requested sequence with lowercase letters indicating the sites of modifications

	Scenario: Export peptide sequence      #where to put this feature
		  Given a user requests the export of a peptide sequence
		  When non-synonymous SNPs exist in the sequence
		  Then export the requested sequence with lowercase letters indicating the sites of modifications and multiple sequences with wild-type and alternate sequences as captured by the SNPs

	Scenario: Viewieng a protein anchored in an experiment
		  Given a user is viewing a protein
		  When they chose that protein from within an experiment
		  Then indicate the sites of modification that are site specific by highlighting the representative box and give them an option to filter out all modifications not captured in that experiment

		 