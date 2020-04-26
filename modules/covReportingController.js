var featuresList1 = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "not (name in ('ORF_1a', 'ORF_1ab')) and featureMetatags.name = 'CODES_AMINO_ACIDS' and featureMetatags.value = true", "name", "displayName", "parent.name"]));

var featuresList2 = glue.tableToObjects(
		glue.command(["list", "feature", "-w", "name in ('ORF_1a', 'ORF_1ab')"]));

var featuresList = featuresList1.concat(featuresList2);


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
	var numSeqs = fastaDocument.nucleotideFasta.sequences.length;
	var processed = 1;
	// apply report generation to each sequence in the set.
	var covReports = _.map(fastaDocument.nucleotideFasta.sequences, function(sequence) {
		glue.setRunningDescription("Sequence analysis for "+processed+"/"+numSeqs+" sequence"+((numSeqs > 1) ? "s" : ""));
		processed++;
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

	
	glue.log("FINE", "covReportingController.initResultMap, result map after placement / lineage assignment", resultMap);
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
		if(featureObj.name == "ORF_1a" || featureObj.name == "ORF_1ab") {
			return; // don't report replacements in these features, rely on NSPs instead.
		}
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
								featureDisplayName: featureObj.displayName,
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
		replacement.displayText = 
			"amino acid replacement in "+replacement.replacement.featureDisplayName+": "+
				refAas.join('/')+replacement.replacement.codonLabel+queryAas.join('/');
		if(knownCovReplacements.length > 0) {
			replacement.knownLink = "replacement/"+knownCovReplacements[0].known_replacement.id;
		}
	});
	return replacementsList;
}

function generateInsertions(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
	var insertionsList = [];
	var insertionsSet = {};
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
				// skip very long insertions 
				if(insertionObj.insertedQryNts.length > 100) {
					return;
				}
				// skip insertions consisting purely of Ns
				if(insertionObj.insertedQryNts.match("^N+$") != null) {
					return;
				}
				insertionObj.featureDisplayName = featureObj.displayName;
				glue.logInfo("insertionObj", insertionObj);
				var insertionKey = insertionObj.refLastNtBeforeIns+":"+insertionObj.insertedQryNts+":"+insertionObj.refFirstNtAfterIns;
				if(insertionsSet[insertionKey] == null) { // don't report same insertion twice for both feature and parent feature.
					insertionsList.push({
						insertion: insertionObj
					});
					insertionsSet[insertionKey] = "yes";
				}
			});
		});
	});
	_.each(insertionsList, function(insertion) {
		if(insertion.insertion.insertionIsCodonAligned) {
			var displayInsertedAas = insertion.insertion.insertedQryAas;
			if(displayInsertedAas.length > 6) {
				displayInsertedAas = displayInsertedAas.substring(0, 3)+"..."+displayInsertedAas.substring(displayInsertedAas.length-3, displayInsertedAas.length);
			}
			insertion.displayText = 
				"codon-aligned insertion in "+insertion.insertion.featureDisplayName+": "+
				insertion.insertion.refLastCodonBeforeIns+"-"+displayInsertedAas+"-"+insertion.insertion.refFirstCodonAfterIns;
		} else {
			var displayInsertedNts = insertion.insertion.insertedQryNts;
			if(displayInsertedNts.length > 6) {
				displayInsertedNts = displayInsertedNts.substring(0, 3)+"..."+displayInsertedNts.substring(displayInsertedNts.length-3, displayInsertedNts.length);
			}
			if(insertion.insertion.insertedQryNts.length % 3 != 0) {
				insertion.displayText = 
					"frameshifting insertion in "+insertion.insertion.featureDisplayName+": ";
			} else {
				insertion.displayText = 
					"non-codon-aligned insertion in "+insertion.insertion.featureDisplayName+": ";
			}
			insertion.displayText += insertion.insertion.refLastNtBeforeIns+"-"+displayInsertedNts+"-"+insertion.insertion.refFirstNtAfterIns;

		}		
		var knownCovInsertions = [];
		var foundCovInsertions = glue.tableToObjects(
				glue.command(["list", "custom-table-row", "cov_nt_insertion", 
					"-w", "variation.featureLoc.feature.name = '"+insertion.insertion.variationFeature+"'"+
					" and last_ref_nt_before = "+insertion.insertion.refLastNtBeforeIns+
					" and first_ref_nt_after = "+insertion.insertion.refFirstNtAfterIns+
					" and inserted_nts = '"+insertion.insertion.insertedQryNts+"'", 
					"id", "display_name", "num_seqs"]));
		_.each(foundCovInsertions, function(foundIns) {
			knownCovInsertions.push({"known_insertion": foundIns});
		});
		insertion.insertion.known_insertions = knownCovInsertions;
		if(knownCovInsertions.length > 0) {
			insertion.knownLink = "insertion/"+knownCovInsertions[0].known_insertion.id;
		}

	});
	return insertionsList;
}

function generateDeletions(queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var comparisonRefName = "REF_MASTER_WUHAN_HU_1";
	var deletionsList = [];
	var deletionsSet = {};
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
				deletionObj.featureDisplayName = featureObj.displayName;
				glue.logInfo("deletionObj", deletionObj);
				var deletionKey = deletionObj.refFirstNtDeleted+":"+deletionObj.refLastNtDeleted;
				if(deletionsSet[deletionKey] == null) { // don't report the same deletion twice for both feature and parent feature.
					deletionsList.push({
						deletion: deletionObj
					});
					deletionsSet[deletionKey] = "yes";
				}
			});
		});
	});
	_.each(deletionsList, function(deletion) {
		if(deletion.deletion.deletionIsCodonAligned) {
			deletion.displayText = 
				"codon-aligned deletion in "+deletion.deletion.featureDisplayName+": ";
			if(deletion.deletion.refFirstCodonDeleted == deletion.deletion.refLastCodonDeleted) {
				deletion.displayText += "codon "+deletion.deletion.refFirstCodonDeleted;
			} else {
				deletion.displayText += "codons "+deletion.deletion.refFirstCodonDeleted+"-"+deletion.deletion.refLastCodonDeleted;
			}
		} else {
			if(deletion.deletion.deletedRefNts.length % 3 != 0) {
				deletion.displayText = 
					"frameshifting deletion in "+deletion.deletion.featureDisplayName+": ";
			} else {
				deletion.displayText = 
					"non-codon-aligned deletion in "+deletion.deletion.featureDisplayName+": ";
			}
			if(deletion.deletion.refFirstNtDeleted == deletion.deletion.refLastNtDeleted) {
				deletion.displayText += "nucleotide "+deletion.deletion.refFirstNtDeleted;
			} else {
				deletion.displayText += "nucleotides "+deletion.deletion.refFirstNtDeleted+"-"+deletion.deletion.refLastNtDeleted;
			}
		}
		var knownCovDeletions = [];
		var foundCovDeletions = glue.tableToObjects(
				glue.command(["list", "custom-table-row", "cov_nt_deletion", 
					"-w", "variation.featureLoc.feature.name = '"+deletion.deletion.variationFeature+"'"+
					" and reference_nt_start = "+deletion.deletion.refFirstNtDeleted+
					" and reference_nt_end = "+deletion.deletion.refLastNtDeleted, 
					"id", "display_name", "num_seqs"]));
			_.each(foundCovDeletions, function(foundDel) {
				knownCovDeletions.push({"known_deletion": foundDel});
			});
		deletion.deletion.known_deletions = knownCovDeletions;
		if(knownCovDeletions.length > 0) {
			deletion.knownLink = "deletion/"+knownCovDeletions[0].known_deletion.id;
		}
	});
	return deletionsList;
}


function generateSingleFastaReport(fastaMap, resultMap, fastaFilePath) {
	
	_.each(_.values(resultMap), function(sequenceResult) {
		if(sequenceResult.isForwardCov) {
			var sequenceId = sequenceResult.id;
			var targetRefName = "REF_MASTER_WUHAN_HU_1";
			var nucleotides = fastaMap[sequenceId].sequence;
			var queryToTargetRefSegs = generateQueryToTargetRefSegs(targetRefName, nucleotides);
			var queryNucleotides = fastaMap[sequenceId].sequence;
			sequenceResult.featuresWithCoverage = generateFeaturesWithCoverage(targetRefName, queryToTargetRefSegs);

			sequenceResult.targetRefName = targetRefName;

			sequenceResult.visualisationHints = visualisationHints(queryNucleotides, targetRefName, queryToTargetRefSegs);

			sequenceResult.primerProbeMismatchReport = 
				primerProbeMismatchReport(fastaFilePath, sequenceId, queryNucleotides, targetRefName, queryToTargetRefSegs);

			var differences = [];
			differences = differences.concat(generateReplacements(queryNucleotides, targetRefName, queryToTargetRefSegs));
			differences = differences.concat(generateInsertions(queryNucleotides, targetRefName, queryToTargetRefSegs));
			differences = differences.concat(generateDeletions(queryNucleotides, targetRefName, queryToTargetRefSegs));

			sequenceResult.differences = differences;
		}
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
		
		glue.setRunningDescription("Lineage assignment for "+numSeqs+" sequence"+((numSeqs > 1) ? "s" : ""));

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
			var resultObj = resultMap[queryLineageResult.queryName];
			resultObj.lineageAssignmentResult = queryLineageResult;
			if(queryLineageResult.placementLineageResults != null) {
				_.each(queryLineageResult.placementLineageResults, function(placementLineageResult) {
					var placement = _.find(resultObj.placements, function(pl) {
						return pl.placementIndex == placementLineageResult.placementIndex;
					});
					if(placementLineageResult.lineages != null && placementLineageResult.lineages.length > 0) {
						placement.bestLineage = placementLineageResult.lineages[placementLineageResult.lineages.length-1];
					}
				});
			}
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

function primerProbeMismatchReport(fastaFilePath, sequenceId, queryNucleotides, targetRefName, queryToTargetRefSegs) {
	var primerProbeMismatchDocument;
	glue.inMode("module/covPrimerProbeMismatch", function() {
		primerProbeMismatchDocument = glue.command({
			"invoke-function": {
				"functionName": "reportSingleFastaFromDocument",
				"document": {
					inputDocument: {
						fastaFilePath: fastaFilePath,
						queryName: sequenceId,
						queryNucleotides: queryNucleotides,
						targetRefName: targetRefName,
						queryToTargetRefSegs: queryToTargetRefSegs
					}
				}
			}
		});
	});
	return primerProbeMismatchDocument;
}
