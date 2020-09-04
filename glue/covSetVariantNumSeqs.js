var processed = 0;
var repIDs = glue.getTableColumn(glue.command(["list", "custom-table-row", "cov_replacement", "id"]), "id");

_.each(repIDs, function(repID) {
	var numSeqs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_replacement_sequence", "-w", 
				"cov_replacement.id = '"+repID+"'"])).length;
	glue.inMode("custom-table-row/cov_replacement/"+repID, function() {
		glue.command(["set", "field", "num_seqs", numSeqs]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set num_seqs for "+processed+" replacements. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});
glue.logInfo("Set num_seqs for "+processed+" replacements. ");
glue.command(["commit"]);
glue.command(["new-context"]);


processed = 0;
var delIDs = glue.getTableColumn(glue.command(["list", "custom-table-row", "cov_deletion", "id"]), "id");

_.each(delIDs, function(delID) {
	var numSeqs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_deletion_sequence", "-w", 
				"cov_deletion.id = '"+delID+"'"])).length;
	glue.inMode("custom-table-row/cov_deletion/"+delID, function() {
		glue.command(["set", "field", "num_seqs", numSeqs]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set num_seqs for "+processed+" deletions. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});
glue.logInfo("Set num_seqs for "+processed+" deletions. ");
glue.command(["commit"]);
glue.command(["new-context"]);

processed = 0;
var ntDelIDs = glue.getTableColumn(glue.command(["list", "custom-table-row", "cov_nt_deletion", "id"]), "id");

_.each(ntDelIDs, function(ntDelID) {
	var numSeqs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_nt_deletion_sequence", "-w", 
				"cov_nt_deletion.id = '"+ntDelID+"'"])).length;
	glue.inMode("custom-table-row/cov_nt_deletion/"+ntDelID, function() {
		glue.command(["set", "field", "num_seqs", numSeqs]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set num_seqs for "+processed+" nt_deletions. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});
glue.logInfo("Set num_seqs for "+processed+" nt_deletions. ");
glue.command(["commit"]);
glue.command(["new-context"]);

processed = 0;
var insIDs = glue.getTableColumn(glue.command(["list", "custom-table-row", "cov_insertion", "id"]), "id");

_.each(insIDs, function(insID) {
	var numSeqs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_insertion_sequence", "-w", 
				"cov_insertion.id = '"+insID+"'"])).length;
	glue.inMode("custom-table-row/cov_insertion/"+insID, function() {
		glue.command(["set", "field", "num_seqs", numSeqs]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set num_seqs for "+processed+" insertions. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});
glue.logInfo("Set num_seqs for "+processed+" insertions. ");
glue.command(["commit"]);
glue.command(["new-context"]);

processed = 0;
var ntInsIDs = glue.getTableColumn(glue.command(["list", "custom-table-row", "cov_nt_insertion", "id"]), "id");

_.each(ntInsIDs, function(ntInsID) {
	var numSeqs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_nt_insertion_sequence", "-w", 
				"cov_nt_insertion.id = '"+ntInsID+"'"])).length;
	glue.inMode("custom-table-row/cov_nt_insertion/"+ntInsID, function() {
		glue.command(["set", "field", "num_seqs", numSeqs]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.logInfo("Set num_seqs for "+processed+" nt_insertions. ");
		glue.command(["commit"]);
		glue.command(["new-context"]);
	}
});
glue.logInfo("Set num_seqs for "+processed+" nt_insertions. ");
glue.command(["commit"]);
glue.command(["new-context"]);



