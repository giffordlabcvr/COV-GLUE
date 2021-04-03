
function ensure_reps(idList) {
	var fields = ["display_name", "codon_label", 
		"codon_label_int", "reference_aa", "reference_nt", "replacement_aa", 
		"radical_hanada_category_i", "radical_hanada_category_ii", "radical_hanada_category_iii", "grantham_distance_double", "grantham_distance_int", "miyata_distance", 
		"parent_feature", "parent_codon_label"];
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
		var featureName = replacementObj["variation.featureLoc.feature.name"];
		// create variation
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

function ensure_nt_deletions(idList) {

	var fields = ["display_name", 
		  "reference_nt_start", "reference_nt_end", "parent_feature"];
	
	_.each(idList, function(ntDelID) {
		var existing = glue.tableToObjects(glue.command(["list", "custom-table-row", 
			"cov_nt_deletion", "-w", "id = '"+ntDelID+"'"]));
		if(existing.length > 0) {
			return; // nt deletion already in DB.
		}
		var jsonCachePath = "build_cache/nt_deletion/"+ntDelID+".json";
		var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
		if(!jsonCacheExists) {
			throw new Error("Sequence build cache JSON refers to nt_deletion that does not have its own build cache JSON");
		}
		var ntDelObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var featureName = ntDelObj["variation.featureLoc.feature.name"];

		// create variation
		var variationName = "cov_nt_del:"+ntDelID;
		glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName, function() {
				glue.command(["create", "variation", variationName, 
				"-t", "nucleotideDeletion", 
				"--nucleotide", ntDelObj.reference_nt_start, ntDelObj.reference_nt_end]);
		});

		// create nt_deletion
		glue.command(["create", "custom-table-row", "cov_nt_deletion", ntDelID]);
		glue.inMode("custom-table-row/cov_nt_deletion/"+ntDelID, function() {
			_.each(fields, function(field) {
				var value = ntDelObj[field];
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

function ensure_deletions(idList) {

	var fields = ["display_name", 
		  "end_codon", "end_codon_int", "parent_end_codon", "parent_feature", "parent_start_codon", 
		  "reference_nt_end", "reference_nt_start", "start_codon", "start_codon_int"];

	_.each(idList, function(delID) {
		var existing = glue.tableToObjects(glue.command(["list", "custom-table-row", 
			"cov_deletion", "-w", "id = '"+delID+"'"]));
		if(existing.length > 0) {
			return; // deletion already in DB.
		}
		var jsonCachePath = "build_cache/deletion/"+delID+".json";
		var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
		if(!jsonCacheExists) {
			throw new Error("Sequence build cache JSON refers to deletion that does not have its own build cache JSON");
		}
		var delObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var featureName = delObj["variation.featureLoc.feature.name"];
		var ntDeletionID = delObj["cov_nt_deletion.id"];

		// create variation
		var variationName = "cov_aa_del:"+delID;
		glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName, function() {
			glue.command(["create", "variation", variationName, 
				"-t", "aminoAcidDeletion", 
				"--labeledCodon", delObj.start_codon, delObj.end_codon]);
		});
		
		// create deletion
		glue.command(["create", "custom-table-row", "cov_deletion", delID]);
		glue.inMode("custom-table-row/cov_deletion/"+delID, function() {
			_.each(fields, function(field) {
				var value = delObj[field];
				if(value != null) {
					glue.command(["set", "field", field, value]);
				};
				glue.command(["set", "link-target", "variation", 
					"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName+
					"/variation/"+variationName]);
				glue.command(["set", "link-target", "cov_nt_deletion", 
					"custom-table-row/cov_nt_deletion/"+ntDeletionID]);
			});
		});

		
	});

}

function ensure_nt_insertions(idList) {

	var fields = ["display_name", "last_ref_nt_before", "first_ref_nt_after", 
	  	"parent_feature", "inserted_nts_length", "inserted_nts"];

	_.each(idList, function(ntInsID) {
		var existing = glue.tableToObjects(glue.command(["list", "custom-table-row", 
			"cov_nt_insertion", "-w", "id = '"+ntInsID+"'"]));
		if(existing.length > 0) {
			return; // nt insertion already in DB.
		}
		var jsonCachePath = "build_cache/nt_insertion/"+ntInsID+".json";
		var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
		if(!jsonCacheExists) {
			throw new Error("Sequence build cache JSON refers to nt_insertion that does not have its own build cache JSON");
		}
		var ntInsObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var featureName = ntInsObj["variation.featureLoc.feature.name"];

		// create variation
		var variationName = "cov_nt_ins:"+ntInsID;
		glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName, function() {
				glue.command(["create", "variation", variationName, 
					"-t", "nucleotideInsertion", 
					"--nucleotide", ntInsObj.last_ref_nt_before, ntInsObj.first_ref_nt_after]);
				glue.inMode("variation/"+variationName, function() {
					glue.command(["set", "metatag", "MIN_INSERTION_LENGTH_NTS", ntInsObj.inserted_nts_length]);
					glue.command(["set", "metatag", "MAX_INSERTION_LENGTH_NTS", ntInsObj.inserted_nts_length]);
				});
		});
		
		// create nt_insertion
		glue.command(["create", "custom-table-row", "cov_nt_insertion", ntInsID]);
		glue.inMode("custom-table-row/cov_nt_insertion/"+ntInsID, function() {
			_.each(fields, function(field) {
				var value = ntInsObj[field];
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

function ensure_insertions(idList) {

	var fields = ["display_name", 
		  "last_codon_before", "last_codon_before_int", "first_codon_after", 
		  "first_codon_after_int", "last_ref_nt_before", "first_ref_nt_after", "parent_feature", 
		  "parent_last_codon_before", "parent_first_codon_after", "inserted_aas_length", "inserted_aas"];

	_.each(idList, function(insID) {
		var existing = glue.tableToObjects(glue.command(["list", "custom-table-row", 
			"cov_insertion", "-w", "id = '"+insID+"'"]));
		if(existing.length > 0) {
			return; // insertion already in DB.
		}
		var jsonCachePath = "build_cache/insertion/"+insID+".json";
		var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
		if(!jsonCacheExists) {
			throw new Error("Sequence build cache JSON refers to insertion that does not have its own build cache JSON");
		}
		var insObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var featureName = insObj["variation.featureLoc.feature.name"];
		var ntInsertionID = insObj["cov_nt_insertion.id"];

		// create variation
		var variationName = "cov_aa_ins:"+insID;
		glue.inMode("reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName, function() {
			glue.command(["create", "variation", variationName, 
			"-t", "aminoAcidInsertion", 
			"--labeledCodon", insObj.last_codon_before, insObj.first_codon_after]);
			glue.inMode("variation/"+variationName, function() {
				glue.command(["set", "metatag", "MIN_INSERTION_LENGTH_AAS", insObj.inserted_aas_length]);
				glue.command(["set", "metatag", "MAX_INSERTION_LENGTH_AAS", insObj.inserted_aas_length]);
			});
		});

		// create insertion
		glue.command(["create", "custom-table-row", "cov_insertion", insID]);
		glue.inMode("custom-table-row/cov_insertion/"+insID, function() {
			_.each(fields, function(field) {
				var value = insObj[field];
				if(value != null) {
					glue.command(["set", "field", field, value]);
				};
				glue.command(["set", "link-target", "variation", 
					"reference/REF_MASTER_WUHAN_HU_1/feature-location/"+featureName+
					"/variation/"+variationName]);
				glue.command(["set", "link-target", "cov_nt_insertion", 
					"custom-table-row/cov_nt_insertion/"+ntInsertionID]);
			});
		});
		
	});

	
}

var whereClause = "analyse_variation = true";


var seqObjs = glue.tableToObjects(glue.command(["list", "sequence", "-w", whereClause, "source.name", "sequenceID"]));

var processed = 0;

glue.logInfo("Populating from JSON build cache files for "+seqObjs.length+" sequences");

_.each(seqObjs, function(seqObj) {
	var sourceName = seqObj["source.name"];
	var sequenceID = seqObj["sequenceID"];

	var jsonCachePath = "build_cache/sequence_results/"+sourceName+"/"+sequenceID+".json";
	var jsonCacheExists = glue.command(["file-util", "file-exists", jsonCachePath]).fileUtilFileExistsResult.exists;
	var cg_reps_from_cache = false;
	var cg_deletions_from_cache = false;
	var cg_insertions_from_cache = false;
	var variation_present = false;

	if(jsonCacheExists) {
		var cacheObj = JSON.parse(glue.command(["file-util", "load-string", jsonCachePath]).fileUtilLoadStringResult.loadedString);
		var nucleotides;
		glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
			nucleotides = glue.command(["show", "nucleotides"]).nucleotidesResult.nucleotides;
		});
		if(nucleotides == cacheObj.nucleotides) {
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
				ensure_nt_deletions(cacheObj.nt_deletion);
				_.each(cacheObj.nt_deletion, function(ntDelID) {
					var linkObjId = ntDelID+":"+sourceName+":"+sequenceID;
					glue.command(["create", "custom-table-row", "cov_nt_deletion_sequence", linkObjId]);
					glue.inMode("custom-table-row/cov_nt_deletion_sequence/"+linkObjId, function() {
						glue.command(["set", "link-target", "cov_nt_deletion", "custom-table-row/cov_nt_deletion/"+ntDelID]);
						glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
					});
				});
				ensure_deletions(cacheObj.deletion);
				_.each(cacheObj.deletion, function(delID) {
					var linkObjId = delID+":"+sourceName+":"+sequenceID;
					glue.command(["create", "custom-table-row", "cov_deletion_sequence", linkObjId]);
					glue.inMode("custom-table-row/cov_deletion_sequence/"+linkObjId, function() {
						glue.command(["set", "link-target", "cov_deletion", "custom-table-row/cov_deletion/"+delID]);
						glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
					});
				});
			}
			if(cacheObj.nt_insertion != null && cacheObj.insertion != null	) {
				cg_insertions_from_cache = true;
				ensure_nt_insertions(cacheObj.nt_insertion);
				_.each(cacheObj.nt_insertion, function(ntInsID) {
					var linkObjId = ntInsID+":"+sourceName+":"+sequenceID;
					glue.command(["create", "custom-table-row", "cov_nt_insertion_sequence", linkObjId]);
					glue.inMode("custom-table-row/cov_nt_insertion_sequence/"+linkObjId, function() {
						glue.command(["set", "link-target", "cov_nt_insertion", "custom-table-row/cov_nt_insertion/"+ntInsID]);
						glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
					});
				});
				ensure_insertions(cacheObj.insertion);
				_.each(cacheObj.insertion, function(insID) {
					var linkObjId = insID+":"+sourceName+":"+sequenceID;
					glue.command(["create", "custom-table-row", "cov_insertion_sequence", linkObjId]);
					glue.inMode("custom-table-row/cov_insertion_sequence/"+linkObjId, function() {
						glue.command(["set", "link-target", "cov_insertion", "custom-table-row/cov_insertion/"+insID]);
						glue.command(["set", "link-target", "sequence", "sequence/"+sourceName+"/"+sequenceID]);
					});
				});
			}
		}
	}
	glue.inMode("sequence/"+sourceName+"/"+sequenceID, function() {
		glue.command(["set", "field", "cg_reps_from_cache", cg_reps_from_cache]);
		glue.command(["set", "field", "cg_deletions_from_cache", cg_deletions_from_cache]);
		glue.command(["set", "field", "cg_insertions_from_cache", cg_insertions_from_cache]);
		glue.command(["set", "field", "variation_present", variation_present]);
	});
	processed++;
	if(processed % 500 == 0) {
		glue.command(["commit"]);
		glue.command(["new-context"]);
		glue.logInfo("Populated sequence-associated data from JSON build cache files for "+processed+" / "+seqObjs.length+" sequences");
	}
});
glue.command(["commit"]);
glue.command(["new-context"]);
glue.logInfo("Populated sequence-associated data from JSON build cache files for "+processed+" / "+seqObjs.length+" sequences");
