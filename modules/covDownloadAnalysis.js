
function downloadAnalysis(document) {
	var covWebReport = document.inputDocument.covWebReport;
	var downloadFileName = document.inputDocument.downloadFileName;
	var downloadFormat = document.inputDocument.downloadFormat;
	var lineFeedStyle = document.inputDocument.lineFeedStyle;

	var tableResult = generateTableResult(covWebReport);
	
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
}

function generateTableResult(covWebReport) {
	glue.logInfo("covWebReport", covWebReport);
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
}