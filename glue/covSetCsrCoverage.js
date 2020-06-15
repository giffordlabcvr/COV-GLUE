var sourceName = 'cov-gisaid';

var memberObjs;
	
glue.inMode("alignment/AL_GISAID_CONSTRAINED", function() {
	memberObjs = glue.tableToObjects(glue.command(["list", "member", "-w", "sequence.source.name = '"+sourceName+"'"]));
});

glue.log("FINEST", "Setting csr_coverage field for "+memberObjs.length+" sequences");

var processed = 0;

_.each(memberObjs, function(memberObj) {
	var sequenceID = memberObj["sequence.sequenceID"];
	var coverage;
	glue.inMode("alignment/AL_GISAID_CONSTRAINED/member/"+memberObj["sequence.source.name"]+"/"+sequenceID, function() {
		var featureCoverageObj = 
			glue.tableToObjects(
					glue.command(["show", "feature-coverage", "-r", "REF_MASTER_WUHAN_HU_1", "-f", "coding_spanning_region", "--excludeNs"])
			)[0];
		coverage = featureCoverageObj.refNtCoveragePct;
	});
	glue.inMode("sequence/"+memberObj["sequence.source.name"]+"/"+sequenceID, function() {
		glue.command(["set", "field", "--noCommit", "csr_coverage", coverage]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.command(["commit"]);
		glue.command(["new-context"]);
		glue.log("FINEST", "Set csr_coverage field for "+processed+" sequences");
	}
});
glue.command(["commit"]);
glue.command(["new-context"]);
glue.log("FINEST", "Set csr_coverage field for "+processed+" sequences");
	