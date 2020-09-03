
function ensure_reps(idList) {
	_.each(idList, function(repID) {
		var existing = glue.tableToObjects(glue.command(["list", "custom-table-row", 
			"cov_replacement", "-w", "id = '"+repID+"'"]));
		if(existing.length > 0) {
			return; // replacement already in DB.
		}
		var jsonCachePath = "build_cache/replacement/"+repID+".json";
		var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
		if(!jsonCacheExists) {
			throw new Error("Sequence build cache JSON refers to replacement that does not have its own build cache JSON");
		}
		var replacementObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		// create variation
		var featureName = replacementObj["variation.featureLoc.feature.name"];
		var variationName = "cov_aa_rpl:"+repID;
		glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName, function() {
			glue.command(["create", "variation", variationName, 
				"-t", "aminoAcidSimplePolymorphism", 
				"--labeledCodon", replacementObj.codon_label, replacementObj.codon_label]);
			glue.inMode("variation/"+variationName, function() {
				glue.command(["set", "metatag", "SIMPLE_AA_PATTERN", replacementObj.replacement_aa]);
				glue.command(["set", "metatag", "MIN_COMBINED_TRIPLET_FRACTION", 0.25]);
			});
		});
		glue.command(["create", "custom-table-row", "cov_replacement", repID]);
		var fields = ["display_name", "codon_label", 
			"codon_label_int", "reference_aa", "reference_nt", "replacement_aa", 
			"radical_hanada_category_i", "radical_hanada_category_ii", "radical_hanada_category_iii", "grantham_distance_double", "grantham_distance_int", "miyata_distance", 
			"parent_feature", "parent_codon_label"];
		glue.inMode("custom-table-row/cov_replacement/"+repID, function() {
			_.each(fields, function(field) {
				var value = replacementObj[field];
				if(value != null) {
					glue.command(["set", "field", field, value]);
				};
				glue.command(["set", "link-target", "variation", 
					"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName+
					"/variation/"+variationName]);
			});
		});
		
		
	});
}


var lineageVersion;

glue.inMode("custom-table-row/cov_project_properties/lineageVersion", function() {
	lineageVersion = glue.command(["show", "property", "value"]).propertyValueResult.value;
});

var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", "analyse_variation = true", "source.name", "sequenceID"]));

var processed = 0;

glue.logInfo("Populating from JSON build cache files for "+seqObjs.length+" sequences");

_.each(seqObjs, function(seqObj) {
	var sourceName = seqObj["source.name"];
	var sequenceID = seqObj["sequenceID"];

	var jsonCachePath = "build_cache/sequence_results/"+sourceName+"/"+sequenceID+".json";
	var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;

	var cg_lineage_from_cache = false;
	var cg_reps_from_cache = false;
	var cg_deletions_from_cache = false;
	var cg_insertions_from_cache = false;

	if(jsonCacheExists) {
		var cacheObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var nucleotides;
		glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
			nucleotides = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
		});
		if(nucleotides == cacheObj.nucleotides) {
			if(cacheObj.lineageVersion == lineageVersion && 
				cacheObj.cov_glue_lineage != null && 
				cacheObj.cov_glue_lw_ratio != null) {

				cg_lineage_from_cache = true;

				glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
						glue.command(["set", "field", "cov_glue_lineage", cacheObj.cov_glue_lineage]);
						glue.command(["set", "field", "cov_glue_lw_ratio", cacheObj.cov_glue_lw_ratio]);
				});
			}
			if(cacheObj.replacement != null) {
				cg_reps_from_cache = true;
				ensure_reps(cacheObj.replacement);
				_.each(cacheObj.replacement, function(repID) {
					var linkObjId = repID+":"+sourceName+":"+sequenceID;
					glue.command(["create", "custom-table-row", "cov_replacement_sequence", linkObjId]);
					glue.inMode("custom-table-row/cov_replacement_sequence/"+linkObjId, function() {
						glue.command(["set", "link-target", "cov_replacement", "custom-table-row/cov_replacement/"+repID]);
						glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
					});
				});
			}
			if(cacheObj.nt_deletion != null && cacheObj.deletion != null) {
				cg_deletions_from_cache = true;
			}
			if(cacheObj.nt_insertion != null && cacheObj.insertion != null	) {
				cg_insertions_from_cache = true;
			}
		}
	}
	glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		glue.command(["set", "field", "cg_lineage_from_cache", cg_lineage_from_cache]);
		glue.command(["set", "field", "cg_reps_from_cache", cg_reps_from_cache]);
		glue.command(["set", "field", "cg_deletions_from_cache", cg_deletions_from_cache]);
		glue.command(["set", "field", "cg_insertions_from_cache", cg_insertions_from_cache]);
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
