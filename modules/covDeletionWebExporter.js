function exportDeletions(downloadFormat, whereClause, sortProperties, fileName, lineFeedStyle) {
	
	var tableResult = previewDeletions(whereClause, sortProperties);
	
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

function previewDeletions(whereClause, sortProperties) {
	var tableResult;
	
	if(sortProperties == null || sortProperties.trim() == "") {
		sortProperties = "id";
	}
	
	tableResult = glue.command(["list", "custom-table-row", "cov_nt_deletion", 
		"-w", whereClause, "-s", sortProperties, 
		"variation.featureLoc.feature.displayName",
		"variation.featureLoc.feature.description",
		"display_name",
		"reference_nt_start",
		"reference_nt_end",
		"cov_deletion.display_name",
		"cov_deletion.start_codon",
		"cov_deletion.end_codon",
		"variation.featureLoc.feature.parent.displayName",
		"cov_deletion.parent_start_codon",
		"cov_deletion.parent_end_codon",
		"num_seqs",
		]);
	
	// rename columns
	tableResult.listResult.column = [
		"genomeRegion",
		"genomeRegionDesc",
		"deletion",
		"refNtStartPosition",
		"refNtEndPosition",
		"aminoAcidDeletion",
		"startCodonNumber",
		"endCodonNumber",
		"altGenomeRegion",
		"altStartCodonNumber",
		"altEndCodonNumber",
		"numSeqs",
	];
	
	return tableResult;
}