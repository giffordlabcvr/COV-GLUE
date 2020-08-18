glue.command(["commit"]);
glue.command(["new-context"]);

// Full set which has not been populated from the cache.
var whereClause = "analyse_variation = true and cg_lineage_from_cache = false";

// test with ~125 sequences
// var whereClause = "sequenceID > 'EPI_ISL_410000' and sequenceID <= 'EPI_ISL_412999'";

// restrict to sequences which have a PANGOLIN lineage, when just interested in 
// comparing the two lineage assignment methods 
// var whereClause = "pang_lineage != null";

glue.command(["multi-unset", "field", "sequence", "-w", whereClause, "cov_glue_lineage"]);
glue.command(["multi-unset", "field", "sequence", "-w", whereClause, "cov_glue_lw_ratio"]);

var numSeqs = glue.command(["count", "sequence", "-w", whereClause]).countResult.count;

var batchSize = 48;
var processed = 0;
var offset = 0;

glue.command(["commit"]);
glue.command(["new-context"]);

glue.log("FINEST", "Assigning lineages for "+numSeqs+" sequences");

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
	glue.log("INFO", "Assigned lineages for "+processed+"/"+numSeqs+" sequences");
}