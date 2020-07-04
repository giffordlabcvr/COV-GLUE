function exportInsertions(downloadFormat, whereClause, sortProperties, fileName, lineFeedStyle) {
	
	var tableResult = previewInsertions(whereClause, sortProperties);
	
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

function previewInsertions(whereClause, sortProperties) {
	var tableResult;
	
	if(sortProperties == null || sortProperties.trim() == "") {
		sortProperties = "id";
	}
	
	tableResult = glue.command(["list", "custom-table-row", "cov_nt_insertion", 
		"-w", whereClause, "-s", sortProperties, 
		"variation.featureLoc.feature.displayName",
		"variation.featureLoc.feature.description",
		"display_name",
		"last_ref_nt_before",
		"first_ref_nt_after",
		"inserted_nts",
		"inserted_nts_length",
		"cov_insertion.display_name",
		"cov_insertion.last_codon_before",
		"cov_insertion.first_codon_after",
		"variation.featureLoc.feature.parent.displayName",
		"cov_insertion.parent_last_codon_before",
		"cov_insertion.parent_first_codon_after",
		"cov_insertion.inserted_aas",
		"cov_insertion.inserted_aas_length",
		"num_seqs",
		]);
	
	// rename columns
	tableResult.listResult.column = [
		"genomeRegion",
		"genomeRegionDesc",
		"insertion",
		"refNtPositionBefore",
		"refNtPositionAfter",
		"insertedNucleotides",
		"insertedNucleotidesLength",
		"aminoAcidInsertion",
		"codonNumberBefore",
		"codonNumberAfter",
		"altGenomeRegion",
		"altCodonNumberBefore",
		"altCodonNumberAfter",
		"insertedAminoAcids",
		"insertedAminoAcidsLength",
		"numSeqs",
	];
	
	return tableResult;
}