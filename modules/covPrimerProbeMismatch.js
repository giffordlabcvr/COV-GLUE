// static map of concrete / ambiguous NT chars to the set of concrete chars which match them.
var ntCharToSubChars = {
    "A": ["A"],
	"C": ["C"],
	"G": ["G"],
	"T": ["T"],
	"R": ["A", "G"],
	"Y": ["C", "T"],
	"K": ["G", "T"],
	"M": ["A", "C"],
	"S": ["C", "G"],
	"W": ["A", "T"],
	"B": ["C", "G", "T"],
	"D": ["A", "G", "T"],
	"H": ["A", "C", "T"],
	"V": ["A", "C", "G"],
	"N": ["A", "C", "G", "T"]
}

// build map concrete / ambiguous NT chars to a list of concrete / ambig chars which match them.
var ntCharToAllowed = {};

_.each(_.pairs(ntCharToSubChars), function(pair) {
	var ntChar = pair[0];
	var subChars = pair[1];
	var allowed;
	if(subChars.length == 1) {
		allowed = [ntChar];
	} else {
		var fullSubChars = subChars.slice(); // copy array
		_.each(_.pairs(ntCharToSubChars), function(pair2) {
			if(_.difference(pair2[1], subChars).length == 0 && fullSubChars.indexOf(pair2[0]) < 0) {
				fullSubChars.push(pair2[0]);
			}
		});
		allowed = fullSubChars;
	}
	ntCharToAllowed[ntChar] = allowed;
});

function reportSingleFasta(fastaFilePath) {
	
	var fastaDoc;
	glue.inMode("module/covFastaUtility", function() {
		fastaDoc = glue.command(["load-nucleotide-fasta", fastaFilePath]);
	});
	if(fastaDoc.nucleotideFasta.sequences.length > 1) {
		throw new Error("FASTA file can have only one sequence");
	}
	if(fastaDoc.nucleotideFasta.sequences.length == 0) {
		throw new Error("FASTA file has no sequences");
	}
	var queryName = fastaDoc.nucleotideFasta.sequences[0].id;
	var queryNucleotides = fastaDoc.nucleotideFasta.sequences[0].sequence;
	
	var resultDoc = singleSequenceReport(fastaFilePath, queryName, queryNucleotides);
	var htmlFilePath = fastaFilePath.substring(0, fastaFilePath.lastIndexOf("."))+"_ppReport.html";

	glue.inMode("module/covPrimerProbeReportTransformer", function() {
		glue.command({"transform-to-file" : {
			commandDocument: resultDoc,
			outputFile: htmlFilePath
		}});
	});
	glue.log("FINEST", "Wrote primer/probe report: "+htmlFilePath);
}

function singleSequenceReport(fastaFilePath, queryName, queryNucleotides) {
	var targetRefName = "REF_MASTER_WUHAN_HU_1";
	var queryToTargetRefSegs = generateQueryToTargetRefSegs(targetRefName, queryNucleotides);
	
	var projectVersion = 
		glue.command(["show","setting","project-version"]).projectShowSettingResult.settingValue;
	var glueVersion = 
		glue.command(["glue-engine","show-version"]).glueEngineShowVersionResult.glueEngineVersion;
	
	var rasVariationMatchDocument;
	
	var matchResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_match:%'");
	var pulloutResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_pullout:%'");
	var insertionResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_insertion:%'");
	var deletionResults = variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, "name like 'cov_pp_seq_deletion:%'");

	var vNameToMatchResult = {};
	_.each(matchResults, function(mr) {vNameToMatchResult[mr.variationName] = mr});
	var vNameToPulloutResult = {};
	_.each(pulloutResults, function(pr) {vNameToPulloutResult[pr.variationName] = pr});
	var vNameToInsertionResult = {};
	_.each(insertionResults, function(ir) {vNameToInsertionResult[ir.variationName] = ir});
	var vNameToDeletionResult = {};
	_.each(deletionResults, function(dr) {vNameToDeletionResult[dr.variationName] = dr});
	
	var assayObjs = glue.tableToObjects(
			glue.command(["list", "custom-table-row", "cov_primer_probe_assay", 
				"-s", "organisation",
				//"-w", "id = 'NIID_2019-nCOV_N'", // Limited test 
				"id", "display_name", "organisation", "url"]));
	
	_.each(assayObjs, function(assayObj) {
		glue.inMode("custom-table-row/cov_primer_probe_assay/"+assayObj.id, function() {
			assayObj.ppObjs = glue.tableToObjects(glue.command(["list", "link-target", "cov_primer_probe", 
				"id", "display_name", "sequence_fwd", "sequence_fwd_regex", "sequence_rev", "sequence_rev_regex", 
				"ref_hits", "sequence_to_scan", "fwd_orientation", "ref_start", "ref_end", "length", 
				"seq_match.name", "seq_pullout.name", "seq_insertion.name", "seq_deletion.name"]));
		});
		_.each(assayObj.ppObjs, function(ppObj) {
			ppObj.seqMatchResult = vNameToMatchResult[ppObj["seq_match.name"]];
			ppObj.seqPulloutResult = vNameToPulloutResult[ppObj["seq_pullout.name"]];
			ppObj.seqInsertionResult = vNameToInsertionResult[ppObj["seq_insertion.name"]];
			ppObj.seqDeletionResult = vNameToDeletionResult[ppObj["seq_deletion.name"]];
		});
		assayObj.anyIssues = false;
		_.each(assayObj.ppObjs, function(ppObj) {
			ppObj.anyIssues = false;
			ppObj.issues = [];
			ppObj.displayIssues = [];
			var insufficientCoverage = [];
			// insufficientCoverageMismatches if alignment doesn't completely cover the region
			// or if there are any Ns in the region.
			if(ppObj.seqMatchResult.sufficientCoverage == false) {
				ppObj.issues.push({ type: "insufficientCoverageMismatches" });
				insufficientCoverage.push("mismatches");
			} else if(ppObj.seqPulloutResult.sufficientCoverage == true &&
					ppObj.seqPulloutResult.matches.length > 0 &&
					ppObj.seqPulloutResult.matches[0].queryNts.indexOf("N") >= 0) {
				ppObj.issues.push({ type: "insufficientCoverageMismatches" });
				insufficientCoverage.push("mismatches");
			}
			if(ppObj.seqInsertionResult.sufficientCoverage == false) {
				ppObj.issues.push({ type: "insufficientCoverageInsertions" });
				insufficientCoverage.push("insertions");
			}
			if(ppObj.seqDeletionResult.sufficientCoverage == false) {
				ppObj.issues.push({ type: "insufficientCoverageDeletions" });
				insufficientCoverage.push("deletions");
			}
			if(insufficientCoverage.length > 0) {
				ppObj.displayIssues.push("Insufficient coverage to determine "+insufficientCoverage.join("/"));
			}
			if(ppObj.seqMatchResult.sufficientCoverage == true && 
					ppObj.seqPulloutResult.matches.length > 0 &&
					ppObj.seqMatchResult.present == false &&
					insufficientCoverage.indexOf("mismatches") < 0) {
				var mismatches = [];
				var pulloutMatch = ppObj.seqPulloutResult.matches[0];
				var pulloutMatchNts = pulloutMatch.queryNts;
				var sequenceToScan = ppObj.sequence_to_scan;
				glue.logInfo("sequenceToScan", sequenceToScan);
				glue.logInfo("pulloutMatchNts", pulloutMatchNts);

				for(var i = 0; i < sequenceToScan.length; i++) {
					var ppSequenceChar = sequenceToScan[i];
					var allowedNts = ntCharToAllowed[ppSequenceChar];
					var queryNt = pulloutMatchNts[i];
					if(allowedNts.indexOf(queryNt) < 0) {
						var refCoord = ppObj.ref_start+i;
						ppObj.issues.push({ type: "mismatch", 
							ppSequenceChar: ppSequenceChar,
							allowedNts: allowedNts,
							queryNt: queryNt,
							refCoord: refCoord,
							queryCoord:pulloutMatch.queryNtStart+i });
						mismatches.push(ppSequenceChar+refCoord+queryNt);
					}
				}
				if(mismatches.length > 0) {
					var displayString = mismatches.length + " " +
						(mismatches.length == 1 ? "mismatch" : "mismatches") + ": " +
							mismatches.join(", ");
					ppObj.displayIssues.push(displayString);
				} else {
					// sanity check.
					glue.logInfo("ntCharToAllowed", ntCharToAllowed);
					glue.logInfo("ppObj.seqMatchResult", ppObj.seqMatchResult);
					glue.logInfo("ppObj.seqPulloutResult", ppObj.seqPulloutResult);
					throw new Error("Could not find mismatched character despite match variation reporting absent.");
				}
			}
			
			if(ppObj.seqDeletionResult.sufficientCoverage == true && ppObj.seqDeletionResult.present == true) {
				var deletions = [];
				_.each(ppObj.seqDeletionResult.matches, function(delMatch) {
					var refFirstNtDeleted = delMatch.refFirstNtDeleted;
					var refLastNtDeleted = delMatch.refLastNtDeleted;
					var deletedRefNts = delMatch.deletedRefNts;
					ppObj.issues.push({ type: "deletion", 
						refFirstNtDeleted: refFirstNtDeleted,
						refLastNtDeleted: refLastNtDeleted,
						qryLastNtBeforeDel: delMatch.qryLastNtBeforeDel,
						qryFirstNtAfterDel: delMatch.qryFirstNtAfterDel,
						deletedRefNts: deletedRefNts });
					var delString;
					if(refFirstNtDeleted == refLastNtDeleted) {
						delString = refFirstNtDeleted;
					} else {
						delString = refFirstNtDeleted+"-"+refLastNtDeleted;
					}
					deletions.push(delString);
				});
				if(deletions.length > 0) {
					var displayString = deletions.length + " " +
						(deletions.length == 1 ? "deletion" : "deletions") + ": " +
							deletions.join(", ");
					ppObj.displayIssues.push(displayString);
				}
			}

			if(ppObj.seqInsertionResult.sufficientCoverage == true && ppObj.seqInsertionResult.present == true) {
				var insertions = [];
				_.each(ppObj.seqInsertionResult.matches, function(insMatch) {
					var refLastNtBeforeIns = insMatch.refLastNtBeforeIns;
					var refFirstNtAfterIns = insMatch.refFirstNtAfterIns;
					var insertedQryNts = insMatch.insertedQryNts;
					ppObj.issues.push({ type: "insertion", 
						refLastNtBeforeIns: refLastNtBeforeIns,
						refFirstNtAfterIns: refFirstNtAfterIns,
						qryFirstInsertedNt: insMatch.qryFirstInsertedNt,
						qryLastInsertedNt: insMatch.qryLastInsertedNt,
						insertedQryNts: insertedQryNts });
					insertions.push(refLastNtBeforeIns+"-"+insertedQryNts+"-"+refFirstNtAfterIns);
				});
				if(insertions.length > 0) {
					var displayString = insertions.length + " " +
						(insertions.length == 1 ? "insertion" : "insertions") + ": " +
							insertions.join(", ");
					ppObj.displayIssues.push(displayString);
				}
			}
			if(ppObj.issues.length > 0) {
				ppObj.anyIssues = true;
				assayObj.anyIssues = true;
			}
		});		
	});
	
	var resultDoc = {
		covPPReport: {
			fastaFilePath: fastaFilePath,
			sequenceID: queryName,
			assays: assayObjs,
			projectVersion: projectVersion,
			glueVersion: glueVersion,
			reportGenerationDate: todaysDate(),
		}
	};
	
	glue.log("FINEST", "resultDoc", resultDoc);
	return resultDoc;
}

function todaysDate() {
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
	return dd + '/' + mm + '/' + yyyy;
}

function variationMatchResults(queryNucleotides, queryToTargetRefSegs, targetRefName, whereClause) {
	var results;
	glue.inMode("module/covFastaSequenceReporter", function() {
		results = glue.command({
			"string-plus-alignment" :{
				"variation":{
					"scan":{
						"fastaString":queryNucleotides,
						"queryToTargetSegs": {
							queryToTargetSegs: {
								alignedSegment: queryToTargetRefSegs
							}
						},
						"whereClause": whereClause,
						"targetRefName":targetRefName,
						"relRefName":targetRefName,
						"linkingAlmtName":"AL_GISAID_UNCONSTRAINED",
						"featureName":"whole_genome",
						"descendentFeatures":"false",
						"excludeAbsent":"false",
						"excludeInsufficientCoverage":"false",
						"showMatchesAsDocument":"true",
						"showMatchesAsTable":"false"
					}
				}
			}
		}).variationScanMatchCommandResult.variations;
	});
	return results;
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
		//glue.log("FINE", "covPrimerProbeMismatch.generateQueryToTargetRefSegs, alignResult", alignResult);
	});
	return alignResult.compoundAlignerResult.sequence[0].alignedSegment;
}


function reportMultiFasta(fastaFilePath) {
	throw new Error("reportMultiFasta is unimplemented");
}

function reportFastaDirectory(fastaDir) {
	throw new Error("reportFastaDirectory is unimplemented");
}