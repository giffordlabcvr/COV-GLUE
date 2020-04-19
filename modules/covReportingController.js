var featuresList = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));


function reportFastaWeb(base64, filePath) {
	glue.log("FINE", "covReportingController.reportFastaWeb invoked");
	var fastaDocument;
	glue.inMode("module/covFastaUtility", function() {
		fastaDocument = glue.command(["base64-to-nucleotide-fasta", base64]);
	});
	var numSequencesInFile = fastaDocument.nucleotideFasta.sequences.length;
	if(numSequencesInFile == 0) {
		throw new Error("No sequences found in FASTA file");
	}
	var maxSequencesWithoutAuth = 50;
	if(numSequencesInFile > maxSequencesWithoutAuth && !glue.hasAuthorisation("covFastaAnalysisLargeSubmissions")) {
		throw new Error("Not authorised to analyse FASTA files with more than "+maxSequencesWithoutAuth+" sequences");
	}
	result = reportDocument({
		reportFastaDocument: {
			"fastaDocument": fastaDocument, 
			"filePath": filePath
		}
	});
	glue.setRunningDescription("Collating report");
	return result;
}

function reportDocument(document) {
	var filePath = document.reportFastaDocument.filePath;
	var fastaDocument = document.reportFastaDocument.fastaDocument;
	var fastaMap = {};
	var resultMap = {};
	var placerResultContainer = {};
	// apply blast recogniser / genotyping together on set, as this is more efficient.
	initResultMap(fastaDocument, fastaMap, resultMap, placerResultContainer);
	// apply report generation to each sequence in the set.
	var covReports = _.map(fastaDocument.nucleotideFasta.sequences, function(sequence) {
		return generateSingleFastaReport(_.pick(fastaMap, sequence.id), _.pick(resultMap, sequence.id), filePath);
	});
	var result = {
		covWebReport:  { 
			results: covReports, 
			placerResult: placerResultContainer.placerResult
		}
	};

	glue.log("FINE", "covReportingController.reportFastaWeb result", result);
	
	return result;
}

/**
 * Entry point for generating a report for a fasta file containing a single sequence.
 */
function reportFasta(fastaFilePath) {
	glue.log("FINE", "covReportingController.reportFasta invoked, input file:"+fastaFilePath);
	// Load fasta and put in a fastaMap
	var fastaDocument;
	glue.inMode("module/covFastaUtility", function() {
		fastaDocument = glue.command(["load-nucleotide-fasta", fastaFilePath]);
	});
	var numSequencesInFile = fastaDocument.nucleotideFasta.sequences.length;
	if(numSequencesInFile == 0) {
		throw new Error("No sequences found in FASTA file");
	}
	if(numSequencesInFile > 1) {
		throw new Error("Please use only one sequence per FASTA file");
	}
	var fastaMap = {};
	var resultMap = {};
	var placerResultContainer = {};
	initResultMap(fastaDocument, fastaMap, resultMap, placerResultContainer);
	var singleFastaReport = generateSingleFastaReport(fastaMap, resultMap, fastaFilePath);
	singleFastaReport.covReport["placerResult"] = placerResultContainer.placerResult;
	return singleFastaReport;
}

function initResultMap(fastaDocument, fastaMap, resultMap, placerResultContainer) {
	glue.log("FINE", "covReportingController.initResultMap fastaDocument:", fastaDocument);
	_.each(fastaDocument.nucleotideFasta.sequences, function(sequenceObj) {
		fastaMap[sequenceObj.id] = sequenceObj;
	});
	// initialise result map.
	var sequenceObjs = _.values(fastaMap);
	_.each(sequenceObjs, function(sequenceObj) {
		resultMap[sequenceObj.id] = { id: sequenceObj.id };
	});
	
	// apply recogniser to fastaMap
	recogniseFasta(fastaMap, resultMap);

	glue.log("FINE", "covReportingController.initResultMap, result map after recogniser", resultMap);

	// apply phylogenetic placement
	placeFasta(fastaMap, resultMap, placerResultContainer);

	glue.log("FINE", "covReportingController.initResultMap, result map after genotyping", resultMap);
}

function generateQueryToTargetRefSegs(targetRefName, nucleotides) {
	var alignerModule;
	glue.inMode("module/covFastaSequenceReporter", function() {
		alignerModule = glue.command(["show", "property", "alignerModuleName"]).moduleShowPropertyResult.propertyValue;
	});
	var alignResult;
	glue.inMode("module/"+alignerModule, function() {
		alignResult = glue.command({align: {
				referenceName: targetRefName,
				sequence: [
				    { 
				    	queryId: "query", 
				    	nucleotides: nucleotides
				    }
				]
			}
		});
		glue.log("FINE", "covReportingController.generateQueryToTargetRefSegs, alignResult", alignResult);
	});
	return alignResult.mafftAlignerResult.sequence[0].alignedSegment;
	
}

function generateFeaturesWithCoverage(targetRefName, queryToTargetRefSegs) {
	var featuresWithCoverage = []; 
	
	_.each(featuresList, function(feature) {
		glue.inMode("module/covFastaSequenceReporter", function() {
			var coveragePercentage = glue.command({
				"alignment-feature-coverage" :{
							"queryToTargetSegs": {
								queryToTargetSegs: {
									alignedSegment: queryToTargetRefSegs
								}
							},
							"targetRefName":targetRefName,
							"relRefName":"REF_MASTER_WUHAN_HU_1",
							"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
							"featureName":feature.name
						}
			}).fastaSequenceAlignmentFeatureCoverageResult.coveragePercentage;
			
			var featureCopy = _.clone(feature);
			featureCopy.coveragePct = coveragePercentage;
			featuresWithCoverage.push(featureCopy);
		});
	});
	return featuresWithCoverage;
}

function generateReplacements(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
	var replacementsList = [];
	_.each(featuresList, function(featureObj) {
		var refAaObjsMap = {};
		glue.inMode("reference/"+comparisonRefName+"/feature-location/"+featureObj.name, function() {
			var refAaObjs = glue.tableToObjects(glue.command(["amino-acid"]));
			_.each(refAaObjs, function(refAaObj) {
				refAaObjsMap[refAaObj.codonLabel] = refAaObj;
			});
		});
		glue.inMode("module/covFastaSequenceReporter", function() {
			var queryAaObjs = glue.tableToObjects(glue.command({
				"string-plus-alignment": { 
					"amino-acid": {
						"fastaString": queryNucleotides,
						"queryToTargetSegs": {
							queryToTargetSegs: {
								alignedSegment: queryToTargetRefSegs
							}
						},
						"targetRefName":targetRefName,
						"relRefName":comparisonRefName,
						"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
						"featureName":featureObj.name
					}
				}
			}));
			_.each(queryAaObjs, function(queryAaObj) {
				// Require no Ns in the codonNts in order to generate a replacement,
				// unless the replacement is unambiguously a single AA residue.
				// This means we are interpreting N as 'unable to sequence' rather than
				// 'equal proportion A, C, G, T' 
				if(queryAaObj.definiteAas != null && queryAaObj.definiteAas != "" && (queryAaObj.definiteAas.length == 1 || queryAaObj.codonNts.indexOf('N') < 0)) {
					refAaObj = refAaObjsMap[queryAaObj.codonLabel];
					if(refAaObj != null && refAaObj.definiteAas != null && refAaObj.definiteAas != "" && 
							refAaObj.definiteAas != queryAaObj.definiteAas) {
						var replacementObj = {
								feature: featureObj.name,
								codonLabel: queryAaObj.codonLabel,
								refAas: refAaObj.definiteAas,
								queryAas: queryAaObj.definiteAas
						};
						glue.logInfo("replacementObj", replacementObj);
						replacementsList.push({
							replacement: replacementObj
						});
					}
				}
			});
		});
	});
	_.each(replacementsList, function(replacement) {
		var refAas = replacement.replacement.refAas.split('');
		var queryAas = replacement.replacement.queryAas.split('');
		var knownCovReplacements = [];
		_.each(refAas, function(refAa) {
			_.each(queryAas, function(queryAa) {
				var knownCovReplacementID = 
					replacement.replacement.feature+":"+refAa+":"+replacement.replacement.codonLabel+":"+queryAa;
				var foundCovReplacements = glue.tableToObjects(
						glue.command(["list", "custom-table-row", "cov_replacement", 
							"-w", "id = '"+knownCovReplacementID+"'", "id", "display_name", "num_seqs"]));
				if(foundCovReplacements.length == 1) {
					knownCovReplacements.push({ "known_replacement": foundCovReplacements[0]});
				}
			});
		});
		replacement.replacement.known_replacements = knownCovReplacements;
	});
	return replacementsList;
}

function generateInsertions(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
	var insertionsList = [];
	_.each(featuresList, function(featureObj) {
		glue.inMode("module/covFastaSequenceReporter", function() {
			var insertionsFound = glue.tableToObjects(glue.command({
				"string-plus-alignment": { 
					"variation": {
						"scan": {
							"fastaString": queryNucleotides,
							"queryToTargetSegs": {
								queryToTargetSegs: {
									alignedSegment: queryToTargetRefSegs
								}
							},
							"targetRefName":targetRefName,
							"relRefName":comparisonRefName,
							"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
							"featureName":featureObj.name,
							"whereClause":"name = 'cov_aa_ins_detect:"+featureObj.name+"'",
							"descendentFeatures": false,
							"excludeAbsent": true,
							"excludeInsufficientCoverage": true,
							"showMatchesAsTable": true,
							"showMatchesAsDocument": false
						}
					}
				}
			}));
			_.each(insertionsFound, function(insertionObj) {
				glue.logInfo("insertionObj", insertionObj);
				insertionsList.push({
					insertion: insertionObj
				});
			});
		});
	});
	_.each(insertionsList, function(insertion) {
		var knownCovInsertions = [];
		if(insertion.insertion.insertionIsCodonAligned) {
			var foundCovInsertions = glue.tableToObjects(
					glue.command(["list", "custom-table-row", "cov_insertion", 
						"-w", "variation.featureLoc.feature.name = '"+insertion.insertion.variationFeature+"'"+
						" and last_codon_before_int = "+insertion.insertion.refLastCodonBeforeIns+
						" and first_codon_after_int = "+insertion.insertion.refFirstCodonAfterIns+
						" and inserted_aas = '"+insertion.insertion.insertedQryAas+"'", 
						"id", "display_name", "num_seqs"]));
				_.each(foundCovInsertions, function(foundIns) {
					knownCovInsertions.push({"known_insertion": foundIns});
				});
		}
		insertion.insertion.known_insertions = knownCovInsertions;
	});
	return insertionsList;
}

function generateDeletions(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
	var deletionsList = [];
	_.each(featuresList, function(featureObj) {
		glue.inMode("module/covFastaSequenceReporter", function() {
			var deletionsFound = glue.tableToObjects(glue.command({
				"string-plus-alignment": { 
					"variation": {
						"scan": {
							"fastaString": queryNucleotides,
							"queryToTargetSegs": {
								queryToTargetSegs: {
									alignedSegment: queryToTargetRefSegs
								}
							},
							"targetRefName":targetRefName,
							"relRefName":comparisonRefName,
							"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
							"featureName":featureObj.name,
							"whereClause":"name = 'cov_aa_del_detect:"+featureObj.name+"'",
							"descendentFeatures": false,
							"excludeAbsent": true,
							"excludeInsufficientCoverage": true,
							"showMatchesAsTable": true,
							"showMatchesAsDocument": false
						}
					}
				}
			}));
			_.each(deletionsFound, function(deletionObj) {
				glue.logInfo("deletionObj", deletionObj);
				deletionsList.push({
					deletion: deletionObj
				});
			});
		});
	});
	_.each(deletionsList, function(deletion) {
		var knownCovDeletions = [];
		if(deletion.deletion.deletionIsCodonAligned) {
			var foundCovDeletions = glue.tableToObjects(
				glue.command(["list", "custom-table-row", "cov_deletion", 
					"-w", "variation.featureLoc.feature.name = '"+deletion.deletion.variationFeature+"'"+
					" and start_codon_int <= "+deletion.deletion.refFirstCodonDeleted+
					" and end_codon_int >= "+deletion.deletion.refLastCodonDeleted, 
					"id", "display_name", "num_seqs"]));
			_.each(foundCovDeletions, function(foundDel) {
				knownCovDeletions.push({"known_deletion": foundDel});
			});
		}
		deletion.deletion.known_deletions = knownCovDeletions;
	});
	return deletionsList;
}


function generateSingleFastaReport(fastaMap, resultMap, fastaFilePath) {
	
	_.each(_.values(resultMap), function(sequenceResult) {
		var targetRefName = "REF_MASTER_WUHAN_HU_1";
		var nucleotides = fastaMap[sequenceResult.id].sequence;
		var queryToTargetRefSegs = generateQueryToTargetRefSegs(targetRefName, nucleotides);
		var queryNucleotides = fastaMap[sequenceResult.id].sequence;
		sequenceResult.featuresWithCoverage = generateFeaturesWithCoverage(targetRefName, queryToTargetRefSegs);

		sequenceResult.targetRefName = targetRefName;

		sequenceResult.visualisationHints = visualisationHints(queryNucleotides, targetRefName, queryToTargetRefSegs);
		
		sequenceResult.replacements = generateReplacements(queryNucleotides, targetRefName, queryToTargetRefSegs);
		sequenceResult.insertions = generateInsertions(queryNucleotides, targetRefName, queryToTargetRefSegs);
		sequenceResult.deletions = generateDeletions(queryNucleotides, targetRefName, queryToTargetRefSegs);
	});
	
	var results = _.values(resultMap);
	
	var covReport = { 
		covReport: {
			sequenceDataFormat: "FASTA",
			filePath: fastaFilePath,
			sequenceResult: results[0]
		}
	};
	addOverview(covReport);

	glue.log("FINE", "covReportingController.generateSingleFastaReport covReport:", covReport);
	return covReport;
}
function visualisationHints(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	// consider the target ref as comparison ref.
	var comparisonReferenceNames = [];
	comparisonReferenceNames.push(targetRefName);
	var seqs = [];
	var comparisonRefs = [];
	
	// eliminate duplicates and enhance with display names.
	_.each(comparisonReferenceNames, function(refName) {
		glue.inMode("reference/"+refName, function() {
			var seqID = glue.command(["show", "sequence"]).showSequenceResult["sequence.sequenceID"];
			if(seqs.indexOf(seqID) < 0) {
				seqs.push(seqID);
				var refDisplayName = glue.command(["show", "property", "displayName"]).propertyValueResult.value;
				if(refDisplayName == null) {
					refDisplayName = "Closest Reference ("+seqID+")";
				}
				comparisonRefs.push({
					"refName": refName,
					"refDisplayName": refDisplayName
				});
			}
		});
	});
	
	var queryDetails = [];
	
	
	return {
		"features": featuresList,
		"comparisonRefs": comparisonRefs,
		"targetReferenceName":targetRefName,
		"queryNucleotides":queryNucleotides,
		"queryToTargetRefSegments": queryToTargetRefSegs,
		"queryDetails": queryDetails
	};
}


/*
 * This function takes a fastaMap of id -> { id, nucleotideFasta }, and a result map of id -> ? 
 * and runs max likelihood placement on the subset of sequences that have been identified as forward nCoV.
 */
function placeFasta(fastaMap, resultMap, placerResultContainer) {
	var placementFastaMap = {};
	_.each(_.values(resultMap), function(resultObj) {
		if(resultObj.isForwardCov && !resultObj.isReverseCov) {
			placementFastaMap[resultObj.id] = fastaMap[resultObj.id];
		} 
	});
	if(!_.isEmpty(placementFastaMap)) {

		var numSeqs = _.values(placementFastaMap).length;
		glue.setRunningDescription("Phylogenetic placement for "+numSeqs+" sequence"+((numSeqs > 1) ? "s" : ""));

		// run the placer and generate a placer result document
		var placerResultDocument;
		glue.inMode("module/covMaxLikelihoodPlacer", function() {
			placerResultDocument = glue.command({
				"place": {
					"fasta-document": {
						"fastaCommandDocument": {
							"nucleotideFasta" : {
								"sequences": _.values(placementFastaMap)
							}
						}
					}
				}
			});
		});
		placerResultContainer.placerResult = placerResultDocument;
		

		
		// list the query summaries within the placer result document
		var placementSummaries;
		glue.inMode("module/covMaxLikelihoodPlacer", function() {
			placementSummaries = glue.tableToObjects(glue.command({
				"list": {
					"query-from-document": {
						"placerResultDocument": placerResultDocument
					}
				}
			}));
		});

		// for each query in the placer results.
		_.each(placementSummaries, function(placementSummaryObj) {
			var queryName = placementSummaryObj.queryName;
			
			var placements;
			
			// list the placements for that query.
			glue.inMode("module/covMaxLikelihoodPlacer", function() {
				placements = glue.tableToObjects(glue.command({
					"list": {
						"placement-from-document": {
							"queryName": queryName,
							"placerResultDocument": placerResultDocument
						}
					}
				}));
			});

			resultMap[queryName].placements = placements;
		});
		
		var lineageAssignmentResultDocument;
		glue.inMode("module/covAssignLineages", function() {
			lineageAssignmentResultDocument = glue.command({
				"invoke-function": {
					"functionName": "assignLineagesFromPlacerDocument",
					"document": placerResultDocument
				}
			});
		});
		_.each(lineageAssignmentResultDocument.covAssignLineagesResult.queryLineageResults, 
				function(queryLineageResult) {
			resultMap[queryLineageResult.queryName].lineageAssignmentResult = queryLineageResult;
		});
	}
}

/*
 * Use the fastaUtility module to reverse complement a FASTA string
 */
function reverseComplement(fastaString) {
	var reverseComplement;
	glue.inMode("module/covFastaUtility", function() {
		var reverseComplementResult = 
			glue.command(["reverse-complement", "string", 
			              "--fastaString", fastaString]);
		reverseComplement = reverseComplementResult.reverseComplementFastaResult.reverseComplement;
	});
	return reverseComplement;
}

/*
 * This function takes a fastaMap of id -> { id, nucleotideFasta }, and a result map of id -> ? 
 * and runs BLAST recogniser, to determine whether the sequence is nCoV, and if so, whether 
 * it is in the forward direction or reverse complement.
 * The result map will have isForwardCov set to true if a forward hit was found, false otherwise
 * It will have isReverseCov set to true if a reverse hit was found, false otherwise
 */
function recogniseFasta(fastaMap, resultMap) {
	var sequenceObjs = _.values(fastaMap);
	_.each(_.values(resultMap), function(resultObj) {
		resultObj.isForwardCov = false;
		resultObj.isReverseCov = false;
	});
	var fastaDocument = {
		"nucleotideFasta" : {
			"sequences" : sequenceObjs
		}
	};
	var numSeqs = sequenceObjs.length;
	glue.setRunningDescription("Sequence recognition for "+numSeqs+" sequence"+((numSeqs > 1) ? "s" : ""));
	var recogniserResults;
	glue.inMode("module/covSequenceRecogniser", function() {
		recogniserResults = glue.tableToObjects(glue.command({
				"recognise": {
					"fasta-document": {
						"fastaCommandDocument": fastaDocument
					}
				}
		}));
	});
	glue.log("FINE", "covReportingController.reportFasta recogniserResults:", recogniserResults);
	_.each(recogniserResults, function(recogniserResult) {
		if(recogniserResult.direction == 'FORWARD') {
			resultMap[recogniserResult.querySequenceId].isForwardCov = true;
		} else if(recogniserResult.direction == 'REVERSE') {
			resultMap[recogniserResult.querySequenceId].isReverseCov = true;
		} 
	});
}

function addOverview(covReport) {
	var today = new Date();
	var dd = today.getDate();
	var mm = today.getMonth()+1; // January is 0!
	var yyyy = today.getFullYear();
	if(dd < 10) {
	    dd = '0'+dd
	} 
	if(mm < 10) {
	    mm = '0'+mm
	} 
	covReport.covReport.reportGenerationDate = dd + '/' + mm + '/' + yyyy;
	covReport.covReport.engineVersion = 
		glue.command(["glue-engine","show-version"]).glueEngineShowVersionResult.glueEngineVersion;
	covReport.covReport.projectVersion = 
		glue.command(["show","setting","project-version"]).projectShowSettingResult.settingValue;
	
}

