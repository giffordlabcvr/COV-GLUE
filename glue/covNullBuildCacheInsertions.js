var dir = "build_cache/sequence_results/cov-gisaid";
var file_names = glue.getTableColumn(glue.command(["file-util", 
	"list-files", "-d", dir]), "fileName");

var processed = 0;
_.each(file_names, function(file_name) {
	if(file_name.endsWith("json")) {
		var cacheObj = JSON.parse(glue.command(["file-util", "load-string", dir+"/"+file_name]).fileUtilLoadStringResult.loadedString);
		delete cacheObj.nt_insertion;
		delete cacheObj.insertion;
		glue.command(["file-util", "save-string", JSON.stringify(cacheObj, null, 2), dir+"/"+file_name]);
	}
	processed++;
	if(processed % 1000 == 0) {
		glue.logInfo("Nullified cached insertions for "+processed+"/"+file_names.length);
	}
});
glue.logInfo("Nullified cached insertions for "+processed+"/"+file_names.length);
