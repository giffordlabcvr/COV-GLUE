glue.command(["multi-unset", "field", "sequence", "-a", "pang_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-a", "pang_representative"]);

// load data from PANGOLIN project lineages table.
var pangSeqObjs;
glue.inMode("module/tabularUtilityCsv", function() {
	pangSeqObjs = glue.tableToObjects(glue.command(["load-tabular", "tabular/lineages.csv"]));
});

_.each(pangSeqObjs, function(pangSeqObj) {
	var matches = pangSeqObj.name.match(/.*(EPI_ISL_\d{6}).*/);
	if(matches == null || matches.length != 2) {
		throw new Error("PANGOLIN lineage data: sequence ID match failed for name value '"+name+"'");
	}
	var seqID = matches[1];
	var numWithSeqId = 
		glue.command(["count", "sequence", 
			"-w", "source.name = 'cov-gisaid' and sequenceID = '"+seqID+"'"]).countResult.count;
	if(numWithSeqId == 0) {
		glue.log("FINEST", "PANGOLIN lineage data: sequence row relates to missing GISAID sequence ID "+seqID);
	} else {
		glue.inMode("sequence/cov-gisaid/"+seqID, function() {
			glue.command(["set", "field", "pang_lineage", pangSeqObj.lineage]);
			glue.command(["set", "field", "pang_representative", 
				pangSeqObj.representative != null && pangSeqObj.representative == "1"]);
		});
	}
	
});