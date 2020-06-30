var sequencesRetrieved = glue.command(["count", "sequence", "-w", "source.name = 'cov-gisaid'"]).countResult.count;
glue.command(["create", "custom-table-row", "cov_project_properties", "sequencesRetrieved"]);
glue.inMode("custom-table-row/cov_project_properties/sequencesRetrieved", function() {
	glue.command(["set", "field", "display_name", "GISAID sequences retrieved"]);
	glue.command(["set", "field", "value", sequencesRetrieved]);
});


var sequencesPassingExclusion = glue.command(["count", "sequence", "-w", "source.name = 'cov-gisaid' and analyse_variation = true"]).countResult.count;
glue.command(["create", "custom-table-row", "cov_project_properties", "sequencesPassingExclusion"]);
glue.inMode("custom-table-row/cov_project_properties/sequencesPassingExclusion", function() {
	glue.command(["set", "field", "display_name", "GISAID sequences which passed exclusion criteria"]);
	glue.command(["set", "field", "value", sequencesPassingExclusion]);
});
