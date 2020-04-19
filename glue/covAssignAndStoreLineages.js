
var whereClause = "true";
// test with ~125 sequences
// var whereClause = "sequenceID > 'EPI_ISL_410000' and sequenceID <= 'EPI_ISL_412999'";

glue.command(["multi-unset", "field", "sequence", "-w", whereClause, "cov_glue_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-w", whereClause, "cov_glue_lw_ratio"]);

var numSeqs = glue.command(["count", "sequence", "-w", whereClause]).countResult.count;

var batchSize = 50;
var processed = 0;
var offset = 0;

glue.log("FINEST", "Assiging lineages for "+numSeqs+" sequences");

while(processed < numSeqs) {
	var batchAssignments;
	glue.inMode("module/covAssignLineages", function() {
		batchAssignments = glue.tableToObjects(glue.command(["invoke-function", 
			"assignLineagesForSequenceBatch", whereClause, offset, batchSize]));
	});
	_.each(batchAssignments, function(batchAssignment) {
		if(batchAssignment.lineage != null) {
			glue.inMode("sequence/"+batchAssignment.queryName, function() {
				glue.command(["set", "field", "--noCommit", "cov_glue_lineage", batchAssignment.lineage]);
				glue.command(["set", "field", "--noCommit", "cov_glue_lw_ratio", batchAssignment.likelihoodWeightRatio]);
			});
		}
	});
	glue.command(["commit"]);
	glue.command(["new-context"]);
	offset += batchSize;
	processed += batchAssignments.length;
	glue.log("FINEST", "Assigned lineages for "+processed+"/"+numSeqs+" sequences");
}