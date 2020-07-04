function exportReplacements(downloadFormat, whereClause, sortProperties, fileName, lineFeedStyle) {
	
	var tableResult = previewReplacements(whereClause, sortProperties);
	
	var saveResult;
	
	if(downloadFormat == "TAB") {
		exportModule = "tabularUtilityTab";
	} else if(downloadFormat == "CSV") {
		exportModule = "tabularUtilityCsv";
	} else {
		throw new Error("Unknown download format '"+downloadFormat+"'");
	}
	
	glue.inMode("module/"+exportModule, function() {
		saveResult = glue.command({
			"save-tabular-web": {
				"lineFeedStyle": lineFeedStyle,
				"fileName": fileName,
				"tabularData": tableResult
			}
		});
	});
	
	return saveResult;
}

function previewReplacements(whereClause, sortProperties) {
	var tableResult;
	
	tableResult = glue.command(["list", "custom-table-row", "cov_replacement", 
		"-w", whereClause, "-s", sortProperties, 
		"variation.featureLoc.feature.displayName",
		"variation.featureLoc.feature.description",
		"display_name",
		"codon_label",
		"variation.featureLoc.feature.parent.displayName",
		"parent_codon_label",
		"reference_nt",
		"reference_aa",
		"replacement_aa",
		"grantham_distance_int",
		"miyata_distance",
		"num_seqs",
		]);
	
	// rename columns
	tableResult.listResult.column = [
		"genomeRegion",
		"genomeRegionDesc",
		"replacement",
		"codonNumber",
		"altGenomeRegion",
		"altCodonNumber",
		"refNtPosition",
		"refAminoAcid",
		"replacementAminoAcid",
		"granthamDistance",
		"miyataDistance",
		"numSeqs",
	];
	
	return tableResult;
}