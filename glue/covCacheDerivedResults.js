

var lineageVersion;

glue.inMode("custom-table-row/cov_project_properties/lineageVersion", function() {
	lineageVersion = glue.command(["show", "property", "value"]).propertyValueResult.value;
});

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", "analyse_variation = true", "source.name", "sequenceID"]));

var processed = 0;

glue.logInfo("Generating JSON build cache files for "+seqObjs.length+" sequences");

_.each(seqObjs, function(seqObj) {
	var sourceName = seqObj["source.name"];
	var sequenceID = seqObj["sequenceID"];
	var cacheObj = {
		lineageVersion: lineageVersion	
	}
	glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		cacheObj.nucleotides = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
		cacheObj.cov_glue_lineage = glue.command(["show", "property", "cov_glue_lineage"]).propertyValueResult.value;
		cacheObj.cov_glue_lw_ratio = glue.command(["show", "property", "cov_glue_lw_ratio"]).propertyValueResult.value;
	});
	var sourceDir = "build_cache/sequence_results/"+sourceName;
	glue.command(["file-util", "make-directories", sourceDir]);
	glue.command(["file-util", "save-string", JSON.stringify(cacheObj, null, 2), sourceDir+"/"+sequenceID+".json"]);
	
	processed++;
	if(processed % 500 == 0) {
		glue.command(["new-context"]);
		glue.logInfo("Generated JSON build cache file for "+processed+" / "+seqObjs.length+" sequences");
	}
});

glue.command(["new-context"]);
glue.logInfo("Generated JSON build cache file for "+processed+" / "+seqObjs.length+" sequences");
