
glue.command(["multi-delete", "cov_project_properties", "-w", "id = 'gisaidTimeStamp'"]);

var timeStampString = glue.command(["file-util", "load-string", "tabular/gisaidTimeStamp.txt"]).fileUtilLoadStringResult.loadedString.trim();

glue.command(["create", "custom-table-row", "cov_project_properties", "gisaidTimeStamp"]);

glue.inMode("custom-table-row/cov_project_properties/gisaidTimeStamp", function() {
	glue.command(["set", "field", "display_name", "Most recent GISAID data update timestamp"]);
	glue.command(["set", "field", "value", timeStampString]);
});
