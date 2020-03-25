glue.command(["multi-unset", "field", "sequence", "-a", "include_in_ref_tree"]);

var clusterRows;

glue.inMode("module/covCdHitEstRunner", function() {
	clusterRows = glue.tableToObjects(glue.command(["generate-clusters", "alignment", "AL_GISAID_CONSTRAINED", 
			"-s", "covReferencePhylogenyColumnsSelector", 
			"-w", "sequence.ref_tree_candidate = true"]));
});

var clusterGroups = _.groupBy(clusterRows, function(cr) {return cr.clusterNumber;});

_.each(_.pairs(clusterGroups), function(pair) {
	var clusterNumber = pair[0];
	var clusterCountries = {};
	_.each(pair[1], function(clusterRow) {
		var sourceName = clusterRow["sequence.source.name"];
		var sequenceID = clusterRow["sequence.sequenceID"];
		glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
			glue.command(["set", "field", "cdhit_cluster", clusterNumber]);
			var country = glue.command(["show", "property", "m49_country.id"]).propertyValueResult.value;
			if(clusterRow.isRepresentative) {
				clusterCountries[country] = "yes";
				glue.logInfo("Including in ref tree cluster "+clusterNumber+
						" CD-HIT representative "+sequenceID+" country "+country);
				glue.command(["set", "field", "include_in_ref_tree", true]);
			} else {
				if(sequenceID == "EPI_ISL_402125") {
					glue.logInfo("Including in ref tree additional sequence for cluster "+
							clusterNumber+" "+sequenceID+ ": country "+country+" (master reference)");
					glue.command(["set", "field", "include_in_ref_tree", true]);
				} else {
					if(clusterCountries[country] == "yes") {
						glue.command(["set", "field", "include_in_ref_tree", false]);
					} else {
						clusterCountries[country] = "yes";
						glue.command(["set", "field", "include_in_ref_tree", true]);
						glue.logInfo("Including in ref tree additional sequence for cluster "+
								clusterNumber+" "+sequenceID+ ": country "+country);
					}
				}
			}
		});
	});
});