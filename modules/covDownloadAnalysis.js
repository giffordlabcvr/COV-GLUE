
function downloadAnalysis(document) {
	var detailLevel = document.inputDocument.detailLevel;
	var covWebReport = document.inputDocument.covWebReport;
	var downloadFileName = document.inputDocument.downloadFileName;
	var downloadFormat = document.inputDocument.downloadFormat;
	var lineFeedStyle = document.inputDocument.lineFeedStyle;

	var tableResult = generateTableResult(detailLevel, covWebReport);
	
	var exportModule;
	if(downloadFormat == "TAB") {
		exportModule = "tabularUtilityTab";
	} else if(downloadFormat == "CSV") {
		exportModule = "tabularUtilityCsv";
	} else {
		throw new Error("Unknown download format '"+downloadFormat+"'");
	}
	var saveResult;
	glue.inMode("module/"+exportModule, function() {
		saveResult = glue.command({
			"save-tabular-web": {
				"lineFeedStyle": lineFeedStyle,
				"fileName": downloadFileName,
				"tabularData": tableResult
			}
		});

	});
	
	return saveResult;
}

function addCovReportRow(covReport, rows, resultType, details) {
	var rowValues = [covReport.filePath, covReport.sequenceResult.id, resultType];
	rows.push({
		"value": rowValues.concat(details)
	});
}

function addResultRowsForCovReport(covReport, rows) {
	var sequenceResult = covReport.sequenceResult;
	var isHCoV19 = sequenceResult.isForwardCov ? "Yes" : "No";
	addCovReportRow(covReport, rows, "isHCoV19", [isHCoV19]);
	var lineage;
	var lwRatio;
	if(sequenceResult.isForwardCov) {
		if(sequenceResult.lineageAssignmentResult != null &&
				sequenceResult.lineageAssignmentResult.bestLineage != null) {
			lineage = sequenceResult.lineageAssignmentResult.bestLineage;
			lwRatio = sequenceResult.lineageAssignmentResult.bestLikelihoodWeightRatio * 100;
		} else {
			lineage = "";
			lwRatio = "";
		}
	} else {
		lineage = "N/A";
		lwRatio = "N/A";
	}
	addCovReportRow(covReport, rows, "lineage", [lineage]);
	addCovReportRow(covReport, rows, "totalLikelihoodWeightRatio", [lwRatio]);
	if(sequenceResult.primerProbeMismatchReport != null) {
		_.each(sequenceResult.primerProbeMismatchReport.covPPReport.assays, function(assay) {
			if(assay.primersWithIssues > 0) {
				_.each(assay.ppObjs, function(ppObj) {
					if(ppObj.numIssues > 0) {
						_.each(ppObj.displayIssues, function(displayIssue) {
							var known = displayIssue.endsWith("*") ? "*KnownIssue" : "";
							addCovReportRow(covReport, rows, "primerProbeIssue", [
								assay.organisation,
								assay.url,
								assay.display_name,
								assay.assay_type,
								ppObj.display_name,
								ppObj.sequence_to_scan,
								ppObj.ref_start+"-"+ppObj.ref_end,
								displayIssue,
								known
							]);
						});
					}
				});
			}
		});
	}
	if(sequenceResult.differences != null) {
		_.each(sequenceResult.differences, function(difference) {
			if(difference.snp != null) {
				addCovReportRow(covReport, rows, "SNP", [
					difference.snp.name
				]);
			} else if(difference.replacement != null) {
				var known = ( difference.replacement.known_replacements != null && 
					difference.replacement.known_replacements.length > 0 ) ? "Known" : "Novel";
				addCovReportRow(covReport, rows, "aaReplacement", [
					known,
					difference.replacement.featureDisplayName,
					difference.replacement.refAas,
					difference.replacement.codonLabel,
					difference.replacement.queryAas
				]);
			} else if(difference.insertion != null) {
				var known = ( difference.insertion.known_insertions != null && 
						difference.insertion.known_insertions.length > 0 ) ? "Known" : "Novel";
				var codonAligned = ( difference.insertion.insertionIsCodonAligned ) ? "Codon-aligned" : "Non-codon-aligned";
					addCovReportRow(covReport, rows, "insertion", [
						known,
						codonAligned,
						difference.insertion.featureDisplayName,
						difference.insertion.refLastNtBeforeIns,
						difference.insertion.insertedQryNts,
						difference.insertion.refFirstNtAfterIns
					]);
			} else if(difference.deletion != null) {
				var known = ( difference.deletion.known_deletions != null && 
						difference.deletion.known_deletions.length > 0 ) ? "Known" : "Novel";
				var codonAligned = ( difference.deletion.deletionIsCodonAligned ) ? "Codon-aligned" : "Non-codon-aligned";
					addCovReportRow(covReport, rows, "deletion", [
						known,
						codonAligned,
						difference.deletion.featureDisplayName,
						difference.deletion.refFirstNtDeleted,
						difference.deletion.refLastNtDeleted
					]);
			}
		});
	}
	
}

function generateTableResult(detailLevel, covWebReport) {
	if(detailLevel == "details") {
		var rows = [];
		_.each(covWebReport.results, function(covReport) {
			addResultRowsForCovReport(covReport.covReport, rows);
		});
		
		return { covDownloadAnalysis : {
			        "column":[
			            "fastaFileName",
			            "sequenceID",
			            "resultType",
			            "detail1",
			            "detail2",
			            "detail3",
			            "detail4",
			            "detail5",
			            "detail6",
			            "detail7",
			            "detail8",
			            "detail9",
			            "detail10",
			        ],
			        "row":rows
		} };
	} else if(detailLevel == "summary") {
		var rows = [];
		_.each(covWebReport.results, function(covReport) {
			rows.push(summaryRowForCovReport(covReport.covReport));
		});
		return { covDownloadAnalysis : {
	        "column":[
	            "fastaFileName",
	            "sequenceID",
	            "isHCoV19",
	            "lineage",
	            "totalLikelihoodWeightRatio",
	            "numDiagnosticsIssues",
	            "numSequencingIssues",
	            "numSnps",
	            "numAaReplacements",
	            "numInsertions",
	            "numDeletions"
	        ],
	        "row":rows
		} };
	}
}


function summaryRowForCovReport(covReport) {
	var sequenceResult = covReport.sequenceResult;
	var isHCoV19 = sequenceResult.isForwardCov ? "Yes" : "No";
	var lineage;
	var lwRatio;
	if(sequenceResult.isForwardCov) {
		if(sequenceResult.lineageAssignmentResult != null &&
				sequenceResult.lineageAssignmentResult.bestLineage != null) {
			lineage = sequenceResult.lineageAssignmentResult.bestLineage;
			lwRatio = sequenceResult.lineageAssignmentResult.bestLikelihoodWeightRatio * 100;
		} else {
			lineage = "";
			lwRatio = "";
		}
	} else {
		lineage = "N/A";
		lwRatio = "N/A";
	}
	var numDiagnosticsIssues = "N/A";
	var numSequencingIssues = "N/A";
	if(sequenceResult.primerProbeMismatchReport != null) {
		numDiagnosticsIssues = sequenceResult.primerProbeMismatchReport.covPPReport.diagnosticsIssues;
		numSequencingIssues = sequenceResult.primerProbeMismatchReport.covPPReport.sequencingIssues;
	}
	var numSnps = "N/A";
	var numAaReplacements = "N/A";
	var numInsertions = "N/A";
	var numDeletions = "N/A";
	
	if(sequenceResult.isForwardCov) {
		numSnps = 0;
		numAaReplacements = 0;
		numInsertions = 0;
		numDeletions = 0;
		if(sequenceResult.differences != null) {
			_.each(sequenceResult.differences, function(difference) {
				if(difference.snp != null) {
					numSnps++;
				} else if(difference.replacement != null) {
					numAaReplacements++;
				} else if(difference.insertion != null) {
					numInsertions++;
				} else if(difference.deletion != null) {
					numDeletions++;
				}
			});
		}
	}
	
	var rowValues = [
		covReport.filePath, 
		covReport.sequenceResult.id,
		isHCoV19,
        lineage,
        lwRatio,
        numDiagnosticsIssues,
        numSequencingIssues,
        numSnps,
        numAaReplacements,
        numInsertions,
        numDeletions];
	return {
		"value": rowValues
	};


	if(sequenceResult.primerProbeMismatchReport != null) {
		_.each(sequenceResult.primerProbeMismatchReport.covPPReport.assays, function(assay) {
			if(assay.primersWithIssues > 0) {
				_.each(assay.ppObjs, function(ppObj) {
					if(ppObj.numIssues > 0) {
						_.each(ppObj.displayIssues, function(displayIssue) {
							var known = displayIssue.endsWith("*") ? "*KnownIssue" : "";
							addCovReportRow(covReport, rows, "primerProbeIssue", [
								assay.organisation,
								assay.url,
								assay.display_name,
								assay.assay_type,
								ppObj.display_name,
								ppObj.sequence_to_scan,
								ppObj.ref_start+"-"+ppObj.ref_end,
								displayIssue,
								known
							]);
						});
					}
				});
			}
		});
	}
	if(sequenceResult.differences != null) {
		_.each(sequenceResult.differences, function(difference) {
			if(difference.snp != null) {
				addCovReportRow(covReport, rows, "SNP", [
					difference.snp.name
				]);
			} else if(difference.replacement != null) {
				var known = ( difference.replacement.known_replacements != null && 
					difference.replacement.known_replacements.length > 0 ) ? "Known" : "Novel";
				addCovReportRow(covReport, rows, "aaReplacement", [
					known,
					difference.replacement.featureDisplayName,
					difference.replacement.refAas,
					difference.replacement.codonLabel,
					difference.replacement.queryAas
				]);
			} else if(difference.insertion != null) {
				var known = ( difference.insertion.known_insertions != null && 
						difference.insertion.known_insertions.length > 0 ) ? "Known" : "Novel";
				var codonAligned = ( difference.insertion.insertionIsCodonAligned ) ? "Codon-aligned" : "Non-codon-aligned";
					addCovReportRow(covReport, rows, "insertion", [
						known,
						codonAligned,
						difference.insertion.featureDisplayName,
						difference.insertion.refLastNtBeforeIns,
						difference.insertion.insertedQryNts,
						difference.insertion.refFirstNtAfterIns
					]);
			} else if(difference.deletion != null) {
				var known = ( difference.deletion.known_deletions != null && 
						difference.deletion.known_deletions.length > 0 ) ? "Known" : "Novel";
				var codonAligned = ( difference.deletion.deletionIsCodonAligned ) ? "Codon-aligned" : "Non-codon-aligned";
					addCovReportRow(covReport, rows, "deletion", [
						known,
						codonAligned,
						difference.deletion.featureDisplayName,
						difference.deletion.refFirstNtDeleted,
						difference.deletion.refLastNtDeleted
					]);
			}
		});
	}
	
}

