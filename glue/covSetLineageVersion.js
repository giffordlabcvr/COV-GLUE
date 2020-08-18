
var lineageVersion = "19th May 2020";

glue.command(["multi-delete", "cov_project_properties", "-w", "id = 'lineageVersion'"]);

glue.command(["create", "custom-table-row", "cov_project_properties", "lineageVersion"]);

glue.inMode("custom-table-row/cov_project_properties/lineageVersion", function() {
	glue.command(["set", "field", "display_name", "Version of the Rambaut, et al. lineage system"]);
	glue.command(["set", "field", "value", lineageVersion]);
});
