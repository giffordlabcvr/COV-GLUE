function exportSequenceMetadata(downloadFormat, whereClause, sortProperties, fileName, lineFeedStyle) {
	
	var tableResult = previewSequenceMetadata(whereClause, sortProperties);
	
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

function previewSequenceMetadata(whereClause, sortProperties) {
	var tableResult;
	
	if(sortProperties == null || sortProperties.trim() == "") {
		sortProperties = "sequenceID";
	}
	
	tableResult = glue.command(["list", "sequence", 
		"-w", whereClause, "-s", sortProperties, 
		"gisaid_virus_name",
		"sequenceID",
		"cov_glue_lineage",
		"gisaid_lineage",
		"gisaid_clade",
		"place_sampled",
		"m49_country.id",
		"m49_country.m49_sub_region.display_name",
		"collection_date",
		]);
	
	// rename columns
	tableResult.listResult.column = [
		"virusName",
		"gisaidID",
		"covGlueLineage",
		"pangolinLineage",
		"gisaidClade",
		"location",
		"countryCode",
		"m49SubRegion",
		"collectionDate",
	];
	
	return tableResult;
}