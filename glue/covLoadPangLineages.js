glue.command(["multi-unset", "field", "sequence", "-a", "pang_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-a", "pang_representative"]);

// load data from PANGOLIN project lineages table.
var pangSeqObjs;
glue.inMode("module/tabularUtilityCsv", function() {
	pangSeqObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/lineages.2020-04-27.csv"]));
});

var processed = 0;

_.each(pangSeqObjs, function(pangSeqObj) {
	var seqID = pangSeqObj["GISAID ID"];
	var representative = pangSeqObj.representative != null && pangSeqObj.representative.trim() == "1";
	var lineage = pangSeqObj.lineage.trim();
	var numWithSeqId = 
		glue.command(["count", "sequence", 
			"-w", "source.name = 'cov-gisaid' and sequenceID = '"+seqID+"'"]).countResult.count;
	if(numWithSeqId == 0) {
		if(representative) {
			throw new Error("Representative "+seqID+" of lineage "+lineage+" no longer exists on GISAID");
		} else {
			glue.log("FINEST", "PANGOLIN lineage data: sequence row relates to missing GISAID sequence ID "+seqID);
		}
	} else {
		glue.inMode("sequence/cov-gisaid/"+seqID, function() {
			glue.command(["set", "field", "pang_lineage", lineage]);
			glue.command(["set", "field", "pang_representative", representative]);
		});
	}
	processed++;
	if(processed % 250 == 0) {
		glue.log("FINEST", "Loaded lineage for "+processed+" sequences");
		glue.command(["new-context"]);
	}
});