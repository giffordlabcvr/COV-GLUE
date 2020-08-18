

var lineageVersion;

glue.inMode("custom-table-row/cov_project_properties/lineageVersion", function() {
	lineageVersion = glue.command(["show", "property", "value"]).propertyValueResult.value;
});

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", "analyse_variation = true", "source.name", "sequenceID"]));

var processed = 0;

glue.logInfo("Populating CoV-GLUE lineage from JSON build cache files for "+seqObjs.length+" sequences");

_.each(seqObjs, function(seqObj) {
	var sourceName = seqObj["source.name"];
	var sequenceID = seqObj["sequenceID"];

	var jsonCachePath = "build_cache/sequence_results/"+sourceName+"/"+sequenceID+".json";
	var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;

	var cg_lineage_from_cache = false;
	if(jsonCacheExists) {
		var cacheObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		if(cacheObj.lineageVersion == lineageVersion) {
			var nucleotides;
			glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
				nucleotides = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
			});
			if(nucleotides == cacheObj.nucleotides) {
				cg_lineage_from_cache = true;
			}
			glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
				if(cacheObj.cov_glue_lineage != null && cacheObj.cov_glue_lw_ratio != null) {
					glue.command(["set", "field", "cov_glue_lineage", cacheObj.cov_glue_lineage]);
					glue.command(["set", "field", "cov_glue_lw_ratio", cacheObj.cov_glue_lw_ratio]);
				}
			});
		}
	}
	glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		glue.command(["set", "field", "cg_lineage_from_cache", cg_lineage_from_cache]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.command(["commit"]);
		glue.command(["new-context"]);
		glue.logInfo("Populated CoV-GLUE lineage from JSON build cache files for "+processed+" / "+seqObjs.length+" sequences");
	}
});

glue.command(["commit"]);
glue.command(["new-context"]);
glue.logInfo("Populated CoV-GLUE lineage from JSON build cache files for "+processed+" / "+seqObjs.length+" sequences");
