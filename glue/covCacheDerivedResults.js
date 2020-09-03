


var variantTypes = [
	{ type: "replacement", 
		  custom_table: "cov_replacement", 
		  fields: ["display_name", "variation.featureLoc.feature.name", "codon_label", 
				"codon_label_int", "reference_aa", "reference_nt", "replacement_aa", 
				"radical_hanada_category_i", "radical_hanada_category_ii", "radical_hanada_category_iii", "grantham_distance_double", "grantham_distance_int", "miyata_distance", 
				"parent_feature", "parent_codon_label"] 
	},

	{ type: "nt_deletion", 
		  custom_table: "cov_nt_deletion", 
		  fields: ["variation.featureLoc.feature.name", "display_name", 
			  "reference_nt_start", "reference_nt_end", "parent_feature"] 
	},

	{ type: "deletion", 
		  custom_table: "cov_deletion", 
		  fields: ["variation.featureLoc.feature.name", "cov_nt_deletion.id", "display_name", 
			  "end_codon", "end_codon_int", "parent_end_codon", "parent_feature", "parent_start_codon", 
			  "reference_nt_end", "reference_nt_start", "start_codon", "start_codon_int" ] 
	},

	{ type: "nt_insertion", 
		  custom_table: "cov_nt_insertion", 
		  fields: ["variation.featureLoc.feature.name", "display_name", "last_ref_nt_before", "first_ref_nt_after", 
			  	"parent_feature", "inserted_nts_length", "inserted_nts"] 
	},

	{ type: "insertion", 
		  custom_table: "cov_insertion", 
		  fields: ["variation.featureLoc.feature.name", "cov_nt_insertion.id", "display_name", 
			  "last_codon_before", "last_codon_before_int", "first_codon_after", 
			  "first_codon_after_int", "last_ref_nt_before", "first_ref_nt_after", "parent_feature", 
			  "parent_last_codon_before", "parent_first_codon_after", "inserted_aas_length", "inserted_aas"] 
	},
	  
]


_.each(variantTypes, function(variantType) {
	var command = ["list", "custom-table-row", variantType.custom_table, "id"];
	command = command.concat(variantType.fields);
	
	var variantObjs = glue.tableToObjects(glue.command(command));

	var variantsDir = "build_cache/"+variantType.type;
	glue.command(["file-util", "make-directories", variantsDir]);

	var processed = 0;

	glue.logInfo("Generating JSON build cache files for "+variantObjs.length+" variants of type "+variantType.type);
	_.each(variantObjs, function(variantObj) {
		glue.command(["file-util", "save-string", JSON.stringify(variantObj, null, 2), variantsDir+"/"+variantObj.id+".json"]);
		processed++;
		if(processed % 500 == 0) {
			glue.logInfo("Generated JSON build cache file for "+processed+" / "+variantObjs.length+" variants of type "+variantType.type);
		}
	});
	glue.logInfo("Generated JSON build cache file for "+processed+" / "+variantObjs.length+" variants of type "+variantType.type);
	glue.command(["new-context"]);
	
});


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
		_.each(variantTypes, function(variantType) {
			cacheObj[variantType.type] = glue.getTableColumn(glue.command(["list", "link-target", variantType.custom_table+"_sequence", variantType.custom_table+".id"]), variantType.custom_table+".id");
		});
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


